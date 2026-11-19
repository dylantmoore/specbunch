# Rebuilding and relinking the LGPL component

The own-code MIT licence permits modification and relinking. CSparse-derived
routines in HiGHS/cuPDLP retain LGPL-2.1-or-later. You may modify those routines,
rebuild the library/plugin, and reverse engineer it to debug those modifications.
No licence term of this package restricts those rights.

Binary wheels, Julia artifacts and Stata net folders include
`corresponding-source.tar.gz` (Stata: `specbunch_source.tar.gz`). It contains our
core, common support, Stata wrapper, build scripts, notices, and the exact pinned
HiGHS/LLVM source archives. No StataCorp SDK source is redistributed.
R source tarballs already contain the complete R wrapper, core and trimmed HiGHS
sources and build rules: edit the sources and run `R CMD INSTALL` to relink.

Extract the archive. With CMake >=3.20, a C++20 compiler, Git and a build tool,
from the extracted specbunch-source directory run (build outside the source):

```
cmake -S bunch/core -B /absolute/build/path -DCMAKE_BUILD_TYPE=Release \
  -DSB_HIGHS_TARBALL="$PWD/third_party/HiGHS-1.8.0.tar.gz" \
  -DSBC_OPENMP_TARBALL="$PWD/third_party/openmp-21.1.8.src.tar.xz" \
  -DSBC_LLVM_CMAKE_TARBALL="$PWD/third_party/cmake-21.1.8.src.tar.xz"
cmake --build /absolute/build/path --config Release --parallel 2
```

To modify HiGHS/CSparse, extract its source archive, make the changes and recreate
the archive; update `SB_HIGHS_SHA256` in the included CMakeLists.txt to the new
archive hash. CMake intentionally verifies source hashes. The same applies to
LLVM archive hashes if modifying that runtime. This verification does not prevent
user modifications. On Windows use a clang-cl developer environment and Ninja.
On macOS set CMAKE_OSX_ARCHITECTURES and CMAKE_OSX_DEPLOYMENT_TARGET as needed.
Use a fresh external build directory to avoid stale objects.

Replace the library inside the installed Python package with the rebuilt
libspecbunch.dylib/libspecbunch.so/specbunch.dll; Julia users can put it in the
artifact lib directory or override the artifact. To rebuild a Stata plugin,
use `-S bunch/stata` and `-DSTATA_PLATFORM=macarm64` (or macintel64, linux64,
win64), retaining the source-archive options above. CMake downloads the unmodified
Stata SDK from StataCorp. Replace the plugin in the net-install directory.
The library ABI and plugin entry points remain available in the included headers.
The system/compiler runtimes and Stata installation are supplied separately.

The Stata source archive additionally contains the two compiled wrapper objects
in prebuilt-stata/. To relink without obtaining or redistributing StataCorp SDK
source, pass both absolute object paths as the semicolon-separated CMake option
`-DSB_STATA_WRAPPER_OBJECTS=/path/plugin.cpp.o;/path/entry.c.o` (quote the whole
argument; Windows uses .obj). This skips SDK download and wrapper compilation,
while rebuilding the modified core and HiGHS. The objects match the architecture
of the accompanying plugin. No SDK source-redistribution permission is assumed.
