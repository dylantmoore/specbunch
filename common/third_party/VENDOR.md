# Private LLVM OpenMP runtime

Pinned LLVM **21.1.8**. Archives are downloaded into the external build tree.

| Archive | SHA-256 |
|---|---|
| openmp-21.1.8.src.tar.xz | `856b023748b41ac7b2c83fd8e9f765ff48a4df2fe6777d2811ef7c7ed8f2f977` |
| cmake-21.1.8.src.tar.xz | `85735f20fd8c81ecb0a09abb0c267018475420e93b65050cc5b7634eab744de9` |

Sources: [OpenMP release archive](https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.8/openmp-21.1.8.src.tar.xz),
[LLVM CMake release archive](https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.8/cmake-21.1.8.src.tar.xz).
The CMake support directory is needed by the standalone OpenMP build. Local
archives or downloads are verified by CMake's URL_HASH before extraction.
No extracted vendor tree or compiled archive is stored in the Drive repository.

`private_openmp.cmake` builds a PIC static `libomp_private` archive in the caller's
external build tree, with private symbol visibility. It is linked separately
into libspecbunch and libspecrd. `SBC_PRIVATE_OPENMP` defaults ON on macOS/Windows,
OFF on Linux. Private builds require the core's static-OpenMP option ON.
Only the package's `sb_*` or `srd_*` C ABI is exported; OpenMP objects and handles
never cross that boundary. No package changes `KMP_DUPLICATE_LIB_OK`.

## Patches

- `libomp-private.patch`, SHA-256
  `9adc53286524ac71d3f3ea05df2fae9ec4ee693f6b5776423783ddbbde8a986b`:
  under `SBC_LIBOMP_PRIVATE`, return from `__kmp_register_library_startup` and
  `__kmp_unregister_library`. Retain startup's x86 system-tick initialization.
  This disables only this copy's duplicate-runtime registration/check/cleanup,
  including LLVM's `__KMP_REGISTERED_LIB_<pid>` environment, shared-memory and
  temp-file mechanism. Without it even hidden static copies abort with Error #15.
  Ordinary OpenMP scheduling, team formation and thread controls are unchanged.
- `libomp-private-windows.patch`, SHA-256
  `d96c6b8cfde45fbe8ba17d54d42d90c2ab2bf8ef3a4b62f897e7ed0e68806de2`:
  allow a private static Windows build through the upstream CMake shared-only
  guard, and skip DLL import-library/export-file generation for static archives.
  Preserve MASM and resource setup; the parent CMake passes the static archive
  itself as the link target expected by upstream Windows micro-test setup.
  MSVC-ABI links also exclude automatic libomp/libiomp5md/vcomp import libraries.
  Compile with `KMP_DYNAMIC_LIB=0` as selected
  by upstream static configuration. This small build-system extension still
  requires execution on Windows; no Windows success is asserted from macOS.

The patches are applied with `git apply --check` then `git apply` to extracted
sources. On a repeated build, a complete `git apply --reverse --check` accepts
an already-applied patch; other mismatches fail. This handles frontend sdist
extraction paths changing while the external build tree is reused. Their combined hash names the build subtree so changed patches cannot
silently reuse old objects. `SBC_LIBOMP_PRIVATE` is defined only for this private
runtime build. OMPT, OMPD, ITT and offloading are disabled. LLVM's Apache 2.0 with
LLVM exceptions license is installed in each wheel as
`<package>/licenses/LLVM-OpenMP-LICENSE.txt` from the release's LICENSE.TXT.

## Linux

Linux builds use GCC's static libgomp, pthreads, static libgcc/libstdc++, and a
version script hiding non-ABI symbols. GCC runtime exception terms apply.
Runtime updates must repeat package coexistence and thread-invariance tests.
