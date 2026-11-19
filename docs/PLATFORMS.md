# Platforms, determinism and threads

Target wheels: macOS arm64 (11.0+) and x86_64 (10.15+), manylinux_2_28 x86_64
and aarch64, and Windows AMD64. R and Julia CI cover Linux, macOS and Windows.
Stata plugins target macarm64, macintel64, linux64 and win64; Linux arm64 has no
Stata target. See REPO_SETUP.md for what has actually been executed locally.

C++20 cores use fixed-order arithmetic, `-ffp-contract=off` and no fast-math.
Within one build, thread changes must preserve numerical bits. Interfaces to
the same core must agree bitwise. Across operating systems, the selected-answer tolerance contract below applies.

Omitted or nonpositive thread counts choose performance cores when detectable,
otherwise physical cores, capped by process availability and readable Linux
cgroup quotas. Positive `OMP_NUM_THREADS` overrides detection. Explicit positive
counts are honored without a core-count cap. `threads_used` records an observed
entry team; some stages remain serial. Builds without OpenMP report one.

macOS and Windows link a patched private static LLVM OpenMP 21.1.8 runtime with
hidden symbols. Linux uses dynamic OpenMP and C++ support libraries; auditwheel bundles libgomp in Linux wheels. Standalone libraries and Stata plugins require the system libgomp.so.1 runtime with GCC.
No package sets `KMP_DUPLICATE_LIB_OK`. R follows its own OpenMP configuration;
its no-runtime fallback requires compiler support for OpenMP SIMD semantics.

## Reproducibility contract

On the same platform and build, results are bitwise identical at any thread count.
Across operating systems, the selected answer agrees on the test fixtures within
relative tolerance 1e-8 (absolute tolerance 1e-12 near zero). Status, refusal reason,
tier, flags, and selected window endpoints must agree exactly; the numerical
comparison covers selected B, sigma, mu (buncher mean), and tau (regression discontinuity),
where applicable. The tolerance is max(1e-12, 1e-8 * max(abs(actual), abs(reference))).

System math libraries differ in the last bit. Diagnostics such as fit p-values and
window ranks can differ slightly, and a sample on a knife edge between two windows
can select a different window. Report the platform with published results.

CI keeps same-platform thread comparisons strictly bitwise. Cross-platform wheel
and Stata-harness jobs fail only for the selected-answer differences above. All
other differences, including candidate-table lengths and unselected windows, are
reported in the job summary and a compressed JSON-lines artifact with every field,
fixture, actual/reference value, and relative difference. Selected-answer fixture
flips are counted explicitly. NPMLE curves have a looser reporting tolerance of
relative 1e-6 or absolute 1e-8 (the larger allowance); exceeding it is reported,
not a cross-platform failure. These tolerances do not change estimator numerics.

The selected-answer tolerance also covers the nonbuncher mean `m_NB` and ETI `saez.e/d` and `continuous.e/d`; rule status, refusal reasons and flags must agree. Other postestimation diagnostics remain reported differences.
