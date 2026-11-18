# Private, source-built LLVM runtime. No runtime symbols may leave the host DLL.
include_guard(GLOBAL)
if(APPLE OR WIN32)
  set(_sbc_private_default ON)
else()
  set(_sbc_private_default OFF)
endif()
option(SBC_PRIVATE_OPENMP "Use a hidden static LLVM OpenMP runtime without process registration" ${_sbc_private_default})
set(SBC_OPENMP_TARBALL "${CMAKE_BINARY_DIR}/downloads/openmp-21.1.8.src.tar.xz" CACHE FILEPATH "Pinned LLVM OpenMP source")
set(SBC_LLVM_CMAKE_TARBALL "${CMAKE_BINARY_DIR}/downloads/cmake-21.1.8.src.tar.xz" CACHE FILEPATH "Pinned LLVM CMake source")
set(_sbc_omp_module "${CMAKE_CURRENT_LIST_DIR}")

function(sbc_link_private_openmp target)
  if(NOT APPLE AND NOT WIN32)
    message(FATAL_ERROR "Private LLVM runtime is supported here only on macOS/Windows; Linux uses static libgomp")
  endif()
  include(ExternalProject)
  find_package(Git REQUIRED)
  if(NOT TARGET sbc_private_omp_build)
    file(SHA256 "${_sbc_omp_module}/libomp-private.patch" patch_hash)
    file(SHA256 "${_sbc_omp_module}/libomp-private-windows.patch" windows_hash)
    string(SHA256 combined_hash "${patch_hash}${windows_hash}")
    string(SUBSTRING "${combined_hash}" 0 12 short_hash)
    set(root "${CMAKE_CURRENT_BINARY_DIR}/private-openmp-${short_hash}")
    foreach(kind IN ITEMS OPENMP LLVM_CMAKE)
      if(EXISTS "${SBC_${kind}_TARBALL}")
        set(${kind}_URL "${SBC_${kind}_TARBALL}")
      elseif(kind STREQUAL "OPENMP")
        set(${kind}_URL "https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.8/openmp-21.1.8.src.tar.xz")
      else()
        set(${kind}_URL "https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.8/cmake-21.1.8.src.tar.xz")
      endif()
    endforeach()
    ExternalProject_Add(sbc_llvm_cmake
      URL "${LLVM_CMAKE_URL}"
      URL_HASH SHA256=85735f20fd8c81ecb0a09abb0c267018475420e93b65050cc5b7634eab744de9
      SOURCE_DIR "${root}/src/cmake" PREFIX "${root}/cmake-stamp"
      DOWNLOAD_EXTRACT_TIMESTAMP ON
      CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "")
    set(extra_args)
    if(WIN32)
      # Upstream Windows micro-test/install setup expects an import-library
      # target. A private static build links the archive itself instead.
      list(APPEND extra_args "-DLIBOMP_IMP_LIB_TARGET=omp"
        "-DLIBOMP_IMP_LIB_FILE=libomp_private${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()
    if(CMAKE_TOOLCHAIN_FILE)
      list(APPEND extra_args "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
    endif()
    if(CMAKE_OSX_SYSROOT)
      list(APPEND extra_args "-DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}")
    endif()
    set(archive "${root}/lib/libomp_private${CMAKE_STATIC_LIBRARY_SUFFIX}")
    ExternalProject_Add(sbc_private_omp_build
      URL "${OPENMP_URL}"
      URL_HASH SHA256=856b023748b41ac7b2c83fd8e9f765ff48a4df2fe6777d2811ef7c7ed8f2f977
      SOURCE_DIR "${root}/src/openmp" BINARY_DIR "${root}/build" PREFIX "${root}/omp-stamp"
      DOWNLOAD_EXTRACT_TIMESTAMP ON
      DEPENDS sbc_llvm_cmake
      PATCH_COMMAND "${CMAKE_COMMAND}" "-DGIT_EXECUTABLE=${GIT_EXECUTABLE}"
        "-DPATCH_DIR=${_sbc_omp_module}" "-DSOURCE_DIR=<SOURCE_DIR>"
        -P "${_sbc_omp_module}/apply_private_patch.cmake"
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DLIBOMP_ENABLE_SHARED=OFF -DSBC_LIBOMP_PRIVATE=ON
        -DLIBOMP_LIB_NAME=libomp_private -DOPENMP_MSVC_NAME_SCHEME=OFF
        -DLIBOMP_INSTALL_ALIASES=OFF -DOPENMP_ENABLE_LIBOMPTARGET=OFF -DOPENMP_ENABLE_OMPT_TOOLS=OFF
        -DLIBOMP_OMPT_SUPPORT=OFF -DLIBOMP_OMPD_SUPPORT=OFF -DLIBOMP_USE_ITT_NOTIFY=OFF
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_VISIBILITY_PRESET=hidden -DCMAKE_VISIBILITY_INLINES_HIDDEN=ON
        "-DCMAKE_CXX_FLAGS=-DSBC_LIBOMP_PRIVATE=1 -ffp-contract=off" "-DCMAKE_C_FLAGS=-DSBC_LIBOMP_PRIVATE=1"
        "-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}"
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        "-DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}" "-DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}"
        "-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=${root}/lib" "-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE=${root}/lib"
        ${extra_args}
      BUILD_COMMAND "${CMAKE_COMMAND}" --build <BINARY_DIR> --target omp --config Release --parallel 4
      INSTALL_COMMAND ""
      BUILD_BYPRODUCTS "${archive}")
    # A normal interface target is visible to callers and inherits dependencies.
    add_library(sbc_private_omp INTERFACE)
    add_dependencies(sbc_private_omp sbc_private_omp_build)
    target_link_libraries(sbc_private_omp INTERFACE "${archive}")
    target_include_directories(sbc_private_omp INTERFACE "${root}/build/runtime/src")
    if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
      target_compile_options(sbc_private_omp INTERFACE -Xpreprocessor -fopenmp)
    elseif(MSVC)
      target_compile_options(sbc_private_omp INTERFACE /clang:-fopenmp)
    else()
      target_compile_options(sbc_private_omp INTERFACE -fopenmp)
    endif()
    if(WIN32)
      target_compile_definitions(sbc_private_omp INTERFACE __KMP_IMP=)
    endif()
    # Notices are shipped inside wheels as well as in sdists.
    install(FILES "${root}/src/openmp/LICENSE.TXT" DESTINATION "${target}/licenses" RENAME LLVM-OpenMP-LICENSE.txt)
  endif()
  target_link_libraries(${target} PRIVATE sbc_private_omp)
  add_dependencies(${target} sbc_private_omp_build)
  set_target_properties(${target} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS OFF)
  if(MSVC)
    # Never satisfy OpenMP references with a driver-selected DLL import library.
    target_link_options(${target} PRIVATE /NODEFAULTLIB:libomp.lib /NODEFAULTLIB:libiomp5md.lib /NODEFAULTLIB:vcomp.lib)
  endif()
  if(MINGW)
    target_link_options(${target} PRIVATE "LINKER:--exclude-all-symbols")
  endif()
endfunction()
