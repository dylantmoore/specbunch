# Joint Voigt numerical checks

The standalone C++20 target `sb_test_joint_voigt` checks nine exact-model transform cases: interior friction, zero Cauchy scale, zero Gaussian variance, and zero/positive/negative jump. Configure `bunch/core` with an external build directory, build this target, then execute it. The fixture CI job runs it before generating the oracle.

`tests/joint_voigt_reference.py --library /external/libspecbunch --output /external/reference.json` independently constructs empirical transforms and covariance weights with NumPy, then compares the selected-window objective and standard error against SciPy constrained least squares. Its four simulated samples cover Gaussian, Voigt, Voigt with a density jump, and Cauchy with a density jump. The objective tolerance is `1e-6 * max(1, reference objective)`; the standard-error tolerance is relative `1e-4` and absolute `1e-8`.

The shared xlang suite checks every named output, including all candidate and transform columns, bitwise within a platform at thread settings 1, 2, 4, 8 and default. Cross-platform comparisons retain the documented selected-answer tolerance contract. All build products, environments, caches, fixtures and validation outputs belong in an external build directory.
