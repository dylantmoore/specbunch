# Platforms, determinism and threads

Target wheels: macOS arm64 (11.0+) and x86_64 (10.15+), manylinux_2_28 x86_64
and aarch64, and Windows AMD64. R and Julia CI cover Linux, macOS and Windows.
Stata plugins target macarm64, macintel64, linux64 and win64; Linux arm64 has no
Stata target. See REPO_SETUP.md for what has actually been executed locally.

C++20 cores use fixed-order arithmetic, `-ffp-contract=off` and no fast-math.
Within one build, thread changes must preserve numerical bits. Interfaces to
the same core must agree bitwise. Across operating systems, math libraries may
change low bits: CI reports bitwise counts and requires relative differences
no larger than 1e-10, with no status, decision, flag or tier changes.

Omitted or nonpositive thread counts choose performance cores when detectable,
otherwise physical cores, capped by process availability and readable Linux
cgroup quotas. Positive `OMP_NUM_THREADS` overrides detection. Explicit positive
counts are honored without a core-count cap. `threads_used` records an observed
entry team; some stages remain serial. Builds without OpenMP report one.

macOS and Windows link a patched private static LLVM OpenMP 21.1.8 runtime with
hidden symbols. Linux uses static libgomp plus static C++ support libraries.
No package sets `KMP_DUPLICATE_LIB_OK`. R follows its own OpenMP configuration;
its no-runtime fallback requires compiler support for OpenMP SIMD semantics.
