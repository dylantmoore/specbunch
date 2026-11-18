# specbunch core: shared library libspecbunch (C ABI in include/specbunch.h). Build OUTSIDE the repo (it is on Drive):
#   cmake -S bunch/core -B ${BUILD_DIR}/bunch -DCMAKE_BUILD_TYPE=Release
#   cmake --build ${BUILD_DIR}/bunch -j
# Options: -DSB_OPENMP=OFF (single-thread build: pragmas ignored, -fopenmp-simd keeps `omp simd` semantics),
#          -DSB_OMP_STATIC=OFF (link the OpenMP runtime dynamically; default on Linux; macOS/Windows default to the private static runtime),
#          -DSBC_PRIVATE_OPENMP=ON (default on macOS/Windows; pinned private LLVM runtime),
#          -DSB_LIBOMP_ROOT=<prefix> (LLVM libomp prefix; set explicitly for non-private macOS builds).
cmake_minimum_required(VERSION 3.20)
include("${CMAKE_CURRENT_LIST_DIR}/../../common/macos_deployment.cmake")
project(specbunch VERSION 0.1.0 LANGUAGES C CXX)

option(SB_OPENMP "Build with OpenMP" ON)
# Distribution libgomp.a uses non-PIC TLS relocations and cannot enter a shared library.
# Linux wheels bundle the dynamic runtime through auditwheel.
if(UNIX AND NOT APPLE)
  set(_omp_static_default OFF)
else()
  set(_omp_static_default ON)
endif()
option(SB_OMP_STATIC "Link the OpenMP runtime statically into libspecbunch" ${_omp_static_default})
set(SB_LIBOMP_ROOT "$ENV{SB_LIBOMP_ROOT}" CACHE PATH "LLVM libomp prefix (macOS)")

if(CMAKE_SOURCE_DIR MATCHES "CloudStorage" AND CMAKE_BINARY_DIR MATCHES "CloudStorage")
  message(FATAL_ERROR "Do not build inside the Drive-synced repo; use -B ${BUILD_DIR}/...")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
# Determinism and portability: no -march=native, no -ffast-math, no FMA contraction; the same flags in every config.
set(CMAKE_CXX_FLAGS_RELEASE "")
# clang-cl accepts GCC optimisation/FP options only through /clang:.
if(MSVC)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(SB_FLAGS /clang:-O3 /DNDEBUG /clang:-ffp-contract=off /clang:-fno-fast-math)
    set(_strict_fp_flag "/clang:-ffp-contract=off")
  else()
    set(SB_FLAGS /O2 /DNDEBUG /fp:strict)
    set(_strict_fp_flag "/fp:strict")
  endif()
else()
  set(SB_FLAGS -O3 -ffp-contract=off -fno-fast-math -fvisibility=hidden -fvisibility-inlines-hidden)
  set(_strict_fp_flag "-ffp-contract=off")
endif()

set(COMMON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../common)
# HiGHS 1.8.0 is the pinned LP solver. Built static, minimal (library target only), PIC, hidden symbols, -ffp-contract=off, outside the repo.
# Source: SB_HIGHS_TARBALL if it exists (default: the copy in pkg_build/third_party), else the GitHub tag tarball;
# either way checked against the SHA-256.
include(ExternalProject)
set(SB_HIGHS_VERSION 1.8.0)
set(SB_HIGHS_SHA256 e184e63101cf19688a02102f58447acc7c021d77eef0d3475ceaceb61f035539)
set(SB_HIGHS_TARBALL "${CMAKE_BINARY_DIR}/downloads/HiGHS-${SB_HIGHS_VERSION}.tar.gz"
    CACHE FILEPATH "Local HiGHS source tarball (used if present)")
if(EXISTS ${SB_HIGHS_TARBALL})
  set(SB_HIGHS_URL ${SB_HIGHS_TARBALL})
else()
  set(SB_HIGHS_URL https://github.com/ERGO-Code/HiGHS/archive/refs/tags/v${SB_HIGHS_VERSION}.tar.gz)
endif()
set(HIGHS_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/highs)
# HiGHS uses the runtime directory for Windows static archives.
if(WIN32)
  set(HIGHS_ARCHIVE_DIR bin)
else()
  set(HIGHS_ARCHIVE_DIR lib)
endif()
set(HIGHS_LIB ${HIGHS_PREFIX}/src/highs_ext-build/${HIGHS_ARCHIVE_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}highs${CMAKE_STATIC_LIBRARY_SUFFIX})
ExternalProject_Add(highs_ext
  URL ${SB_HIGHS_URL}
  URL_HASH SHA256=${SB_HIGHS_SHA256}
  PREFIX ${HIGHS_PREFIX}
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DFAST_BUILD=ON
             -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_BINDIR=bin
             -DZLIB=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_VISIBILITY_PRESET=hidden
             -DCMAKE_VISIBILITY_INLINES_HIDDEN=ON "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} ${_strict_fp_flag}"
             -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
             -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
             -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
             -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
             -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
  BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target highs --config Release --parallel 4
  INSTALL_COMMAND ""
  BUILD_BYPRODUCTS ${HIGHS_LIB}
  DOWNLOAD_EXTRACT_TIMESTAMP ON)
set(HIGHS_INC ${HIGHS_PREFIX}/src/highs_ext/src ${HIGHS_PREFIX}/src/highs_ext-build)

add_library(specbunch SHARED src/abi_common.cpp src/abi_mass_gaussian.cpp src/abi_mass_gaussian_jump.cpp src/abi_bmean.cpp src/npmle.cpp vendor/c2/kernel_puf.cpp vendor/c2/probit.cpp)
add_dependencies(specbunch highs_ext)
target_include_directories(specbunch BEFORE PRIVATE ${COMMON_DIR}/shim/portable)
target_include_directories(specbunch PRIVATE ${HIGHS_INC})
target_link_libraries(specbunch PRIVATE ${HIGHS_LIB})
target_include_directories(specbunch PRIVATE include ${COMMON_DIR} src vendor)
target_compile_definitions(specbunch PRIVATE SBC_BUILDING)
target_compile_options(specbunch PRIVATE ${SB_FLAGS} -Wno-unused-function -Wno-unused-variable)
set_target_properties(specbunch PROPERTIES OUTPUT_NAME specbunch)
# Export only the C ABI (sb_*): keeps a statically linked OpenMP runtime's symbols private to the library.
if(APPLE)
  target_link_options(specbunch PRIVATE "LINKER:-exported_symbol,_sb_*" "LINKER:-dead_strip")
elseif(UNIX)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gen/specbunch.map "{ global: sb_*; local: *; };\n")
  target_link_options(specbunch PRIVATE "LINKER:--version-script=${CMAKE_CURRENT_BINARY_DIR}/gen/specbunch.map" "LINKER:--gc-sections")
  target_compile_options(specbunch PRIVATE -ffunction-sections -fdata-sections)
endif()

include(${COMMON_DIR}/third_party/private_openmp.cmake)
set(SB_OMP_DESC "off")
if(SB_OPENMP AND SBC_PRIVATE_OPENMP AND (APPLE OR WIN32))
  if(NOT SB_OMP_STATIC)
    message(FATAL_ERROR "SBC_PRIVATE_OPENMP requires SB_OMP_STATIC=ON")
  endif()
  sbc_link_private_openmp(specbunch)
  set(SB_OMP_DESC "LLVM 21.1.8 private static; registry disabled")
elseif(SB_OPENMP)
  if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND EXISTS ${SB_LIBOMP_ROOT}/include/omp.h)
    target_compile_options(specbunch PRIVATE -Xpreprocessor -fopenmp)
    target_include_directories(specbunch PRIVATE ${SB_LIBOMP_ROOT}/include)
    if(SB_OMP_STATIC AND EXISTS ${SB_LIBOMP_ROOT}/lib/libomp.a)
      target_link_libraries(specbunch PRIVATE ${SB_LIBOMP_ROOT}/lib/libomp.a)
      set(SB_OMP_DESC "libomp static")
    else()
      target_link_libraries(specbunch PRIVATE ${SB_LIBOMP_ROOT}/lib/libomp.dylib)
      set(SB_OMP_DESC "libomp dynamic")
    endif()
  else()
    find_package(Threads REQUIRED)
    find_package(OpenMP COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
      if(SB_OMP_STATIC AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libgomp.a OUTPUT_VARIABLE GOMP OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT EXISTS "${GOMP}")
          message(FATAL_ERROR "Static libgomp.a is required")
        endif()
        separate_arguments(omp_flags NATIVE_COMMAND "${OpenMP_CXX_FLAGS}")
        target_compile_options(specbunch PRIVATE ${omp_flags})
        target_link_libraries(specbunch PRIVATE "${GOMP}" Threads::Threads)
        target_link_options(specbunch PRIVATE -static-libgcc -static-libstdc++)
        set(SB_OMP_DESC "libgomp static")
      else()
        if(UNIX AND NOT APPLE AND SB_OMP_STATIC)
          message(FATAL_ERROR "Linux static wheel builds require GCC and libgomp.a")
        endif()
        target_link_libraries(specbunch PRIVATE OpenMP::OpenMP_CXX)
        set(SB_OMP_DESC "OpenMP::OpenMP_CXX")
      endif()
    else()
      set(SB_OPENMP OFF)
      message(WARNING "OpenMP not found: single-thread build")
    endif()
  endif()
endif()
if(NOT SB_OPENMP)
  target_include_directories(specbunch BEFORE PRIVATE ${COMMON_DIR}/shim/noomp)
  target_compile_options(specbunch PRIVATE -fopenmp-simd)
  set(SB_OMP_DESC "off (-fopenmp-simd)")
endif()
string(REPLACE ";" " " SB_FLAGS_STR "${SB_FLAGS}")
set(SB_VERSION_STRING "${PROJECT_VERSION}.dev1")
set(SB_BUILD_FLAGS "-std=c++20 ${SB_FLAGS_STR}; omp=${SB_OMP_DESC}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gen/sb_build_config.h.in "#define SB_VERSION_STRING \"@SB_VERSION_STRING@\"\n#define SB_BUILD_FLAGS \"@SB_BUILD_FLAGS@\"\n")
configure_file(${CMAKE_CURRENT_BINARY_DIR}/gen/sb_build_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/gen/sb_build_config.h @ONLY)
target_include_directories(specbunch PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/gen)
target_compile_definitions(specbunch PRIVATE SB_HAVE_BUILD_CONFIG)
message(STATUS "specbunch: flags='${SB_FLAGS_STR}' openmp=${SB_OMP_DESC}")
