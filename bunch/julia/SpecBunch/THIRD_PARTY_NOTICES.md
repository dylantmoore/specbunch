# Third-party notices for specbunch

Dylan Moore's own code, including the research estimator files under `vendor/`,
is licensed under MIT (root LICENSE). This does not relicense third-party code.
Upstream copyright and licence headers remain intact.

| Component | Licence | Where present | Compliance material |
|---|---|---|---|
| HiGHS 1.8.0 | MIT | Statically linked core in wheels, Julia libraries, Stata plugins; trimmed source in R tarballs | HiGHS licence and upstream headers retained. |
| filereaderlp | MIT | HiGHS LP reader / R source | Michael Feldmeier copyright and full licence. |
| pdqsort | zlib | HiGHS sorting / R source | Orson Peters copyright and exact zlib licence (not MIT). |
| zstr 1.0.5 | MIT | R source only; compressed input disabled (ZLIB=OFF) | Matei David / Ontario Institute for Cancer Research licence. zlib itself is not linked. |
| cuPDLP-C 0.3.0 | MIT, except CSparse-derived files | HiGHS CPU solver / R source; no CUDA | COPT-Public licence; HiGHS's adaptation notes retained. |
| CSparse-derived routines (cuPDLP cupdlp_cs.c/.h; IPX sparse_utils; BASICLU lu_dfs) | LGPL-2.1-or-later | HiGHS / R source | Timothy A. Davis notice, LGPL text and provenance; binary redistribution requires corresponding source and relinking materials. |
| IPX | MIT option selected | HiGHS interior-point solver / R source | IPX, Copyright (c) 2018-2020 ERGO-Code. Used in specbunch under the MIT license. |
| BASICLU | MIT option selected | HiGHS sparse factorization / R source | BASICLU, Copyright (c) 2016-2020 ERGO-Code. Used in specbunch under the MIT license. |
| Catch2 2.13.8, embedded Clara 1.1.5 and TextFlowCpp | BSL-1.0 | Downloaded HiGHS source archive only; not compiled or included in trimmed R source | Two Blue Cubes Ltd / Catch2 authors / Phil Nash attribution; Boost licence and embedded headers retained with source. |
| LLVM OpenMP 21.1.8 | Apache-2.0 WITH LLVM-exception | Private static runtime in macOS/Windows wheels, Julia libraries and Stata plugins; not the R source package runtime | Full upstream licence, NOTICE and both modification patches. LLVM CMake support is build-only under the same terms. |
| GCC libgomp (toolchain version) | GPL-3.0-or-later WITH GCC-exception-3.1 | Linux wheels copy the shared runtime using auditwheel; Linux standalone libraries and R use the installed runtime | GPL and runtime-exception texts; see runtime source distribution requirements below. |
| Stata Plugin Interface 3.0 | StataCorp terms, not MIT | Compiled into Stata plugins only; SDK source fetched at build time | Preserve copyright; no SDK source redistributed. See finding below. |

## Sources and scope

Pinned build sources and hashes are recorded in the core CMakeLists.txt and
common/third_party/VENDOR.md. HiGHS components above refer to **1.8.0**, not the
current HiGHS release. specrd does not use HiGHS or any of its components.
Licence texts are in `licenses/`; LLVM modifications are also recorded in NOTICE.

## StataCorp finding

The official https://www.stata.com/plugins/ documentation expressly describes
distributing compiled plugins while keeping their source private. The downloaded
https://www.stata.com/plugins/stplugin.c and
https://www.stata.com/plugins/stplugin.h contain StataCorp LP copyright notices
but no express source-redistribution licence. We therefore do not redistribute
those SDK source files: CMake fetches the unmodified, SHA-256-pinned files from
StataCorp at build time. Plugins retain the acknowledgement: stplugin.c copyright
(c) 2003, 2006, 2015 StataCorp LP; stplugin.h copyright (c) 2003, 2004, 2006, 2015
StataCorp LP. This is not a claim that StataCorp's SDK is MIT licensed.

## Runtime source distribution requirements

The GCC runtime exception permits eligible compiled application code under its
own licence; it does not remove the GPL requirements for redistributing libgomp
itself. A binary release that bundles libgomp must provide access to the matching
corresponding source (including distributor patches), alongside the licence and
exception. Linux wheel builds identify the exact runtime by RPM ownership and SHA-256,
and verify access to its matching source RPM before building. The resulting
GCC-CORRESPONDING-SOURCE.txt is included in wheel licence metadata. Its URL
provides GCC source, distributor patches and the build specification at no charge.
Retain equivalent source access when mirroring or redistributing these wheels.
System-provided runtimes are not redistributed by this repository.

## Dependencies not bundled

NumPy and optional matplotlib are installed separately by Python package managers;
R base, graphics and stats are supplied by R; Julia and optional Plots.jl are
installed separately. Stata itself, operating-system C/C++ libraries, Windows
runtime DLLs, compiler tools, CMake, Ninja, scikit-build-core, pytest and Julia's
package manager are not redistributed as package payloads. Their separately
installed distributions retain their own terms. Compiler runtime code incorporated
into binaries remains governed by the applicable runtime licence/exception.

## LGPL corresponding source

Binary distributions include corresponding-source.tar.gz (Stata:
specbunch_source.tar.gz), with complete library/plugin source, pinned dependency
source archives and rebuild/relink instructions in docs/RELINKING.md. R tarballs
already include the complete corresponding R/core/HiGHS sources. Users may modify
and relink the LGPL portion, and reverse engineer the binary for debugging such
modifications. The own-code MIT licence does not restrict these rights.

## LLVM source-only components

The complete OpenMP source archive also contains Intel ITT Notify under
BSD-3-Clause (licenses/Intel-ITT.txt); it is disabled and not linked. The full
LLVM-OpenMP.txt retains the legacy Intel MIT and University of Illinois/NCSA
notices as well as Apache-2.0 WITH LLVM-exception and the Intel grant.

## Compiler-supplied code

Template instantiations and compiler support routines can be incorporated into
all native artifacts. GCC libstdc++/libgcc use GPL-3.0-or-later WITH
GCC-exception-3.1; LLVM libc++/compiler-rt and the Microsoft C++ Standard Library
use Apache-2.0 WITH LLVM-exception. Their exceptions permit eligible compiled
application code under its own licence. The supplied licence texts retain the
legacy LLVM MIT/NCSA notices too. Microsoft STL attribution: Copyright Microsoft
Corporation; LLVM attribution: LLVM Project contributors. These are compiler
inputs, not additional vendored project source trees. The dynamically imported
OS C/C++ and Windows VC/UCRT runtime libraries themselves are not copied into our
artifacts. GCC link inputs and their source RPMs are recorded separately from
the actually bundled ELF libgomp runtime in Linux wheel metadata.

Upstream texts: https://github.com/microsoft/STL/blob/main/LICENSE.txt,
https://github.com/llvm/llvm-project/blob/llvmorg-21.1.8/libcxx/LICENSE.TXT,
https://github.com/llvm/llvm-project/blob/llvmorg-21.1.8/compiler-rt/LICENSE.TXT.
