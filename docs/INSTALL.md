# Installation details

The copy-paste routes for users are in the [README](../README.md#installation). This page covers what they need, what to do on Windows, common failures, and the developer builds used by CI.

Nothing is on PyPI, CRAN, the Julia General registry or SSC yet; every route compiles the C++20 core from source. Prebuilt binaries (Python wheels, Stata plugins) are planned, so that no compiler will be needed.

## What the build needs

- A C++20 compiler: Apple clang from the Xcode Command Line Tools on macOS, GCC or Clang on Linux, `clang-cl` (LLVM) with the Visual Studio Build Tools on Windows (Rtools for R).
- CMake 3.20 or newer for the Julia and Stata routes. The Python build gets CMake and Ninja from PyPI if they are missing; the R package needs no CMake.
- Git, and an internet connection during the build: CMake downloads HiGHS 1.8.0 (the LP solver) and, on macOS and Windows, the LLVM OpenMP 21.1.8 runtime, both from GitHub and checked against pinned SHA-256 hashes. The R package bundles its sources and downloads nothing.
- Build directories must not be inside a cloud-synced folder (Google Drive, Dropbox, OneDrive); the Python build refuses such paths.

## Threads and OpenMP

The Python wheel, the Julia library and the Stata plugin on macOS and Windows link a private static copy of the LLVM OpenMP runtime, so nothing extra is installed. On Linux they use the system `libgomp` (shipped with GCC). R uses R's own OpenMP configuration: CRAN's macOS R has none, so the R package runs single-threaded there. Results are bitwise identical at any thread count on the same platform, so a single-threaded build gives the same numbers, only more slowly.

R compiles about 180 files one at a time by default. To use several cores, run `Sys.setenv(MAKEFLAGS = "-j4")` before `remotes::install_github(...)`.

## Windows

Windows is the least tested platform, and the steps below follow the CI configuration; they have not been run by hand.

- **R:** install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R version, then use the README command.
- **Python, Julia, Stata:** install the Visual Studio 2022 Build Tools (workload "Desktop development with C++"), LLVM (which provides `clang-cl`; `winget install LLVM.LLVM`), CMake, Ninja and Git. Open the "x64 Native Tools Command Prompt for VS 2022", make sure `C:\Program Files\LLVM\bin` is on `PATH`, and set the compiler and generator the CI uses:

  ```bat
  set CC=clang-cl
  set CXX=clang-cl
  set CMAKE_GENERATOR=Ninja
  ```

  Then run the README commands from that prompt. For Julia, the core library is `specbunch-build\specbunch.dll`; for Stata use `--platform win64 --build-root "%CD%\specbunch-stata"`.

## Troubleshooting

- **The build cannot download HiGHS or OpenMP.** Check access to `github.com`. Offline builds can pass local copies of the same tarballs with `-DSB_HIGHS_TARBALL=...`, `-DSBC_OPENMP_TARBALL=...` and `-DSBC_LLVM_CMAKE_TARBALL=...` (the hashes are still checked).
- **Stata build fails with "unable to make temporary file".** Pass an absolute `--build-root`.
- **Stata cannot load the plugin on macOS.** A locally built plugin loads normally. If you copied a `.plugin` file that was downloaded through a browser, remove the quarantine flag: `xattr -d com.apple.quarantine specbunch_macarm64.plugin` in the PLUS `s/` folder.
- **Julia: `core library not found`.** Run the snippet from the directory that contains `specbunch` and `specbunch-build`, after the CMake build has finished.
- **Python: matplotlib warning when fitting.** Install the plot extra (`specbunch[plot]` in the README command), or pass `plot=False`.

Please report anything else at <https://github.com/dylantmoore/specbunch/issues>.

## Developer and CI builds

These are the commands CI runs; all outputs go to the `--root` directory, outside the checkout.

- Python: `python bunch/python/tools/build.py --root /path/to/build`, then install the wheel from `/path/to/build/py/dist`.
- R: run `R CMD build /path/to/repo/bunch/r` from an external build directory, then `R CMD INSTALL` the resulting archive.
- Julia: `python tools/ci.py julia --root /path/to/build` (builds the core, then tests the package in an isolated depot).
- Stata: `python tools/ci.py stata --root /path/to/build --platform macarm64` (choose the platform matching Stata; this tests through the no-Stata harness and leaves a `net install` folder in `/path/to/build/stata/net/specbunch`).
