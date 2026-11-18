# Spline validation

All inputs are simulated. `reference.py` independently constructs SciPy B-splines, analytic Gaussian integrals, analytic Cauchy spline integrals with high-order Gaussian expectation, constrained optimization, and a Householder-QR EDF calculation. It does not call the package or its C++ core.

`validate.py --output /external/build/validation.json` compares the default estimator on eight cases: 50,000 and 1,000,000 observations, Gaussian and Voigt friction, each with and without a latent density jump. It requires the same selected lambda (relative tolerance 1e-10), absolute log-likelihood error below .003, EDF error below 1e-4, maximum histogram-column error below 1e-6, and maximum curve-column error below 2e-6. Every returned numerical bit is compared across 1, 2, 4 and 8 threads. `--case` runs one numbered case; `--exclude-case` supports disjoint validation batches.

`test_grid.py --library /external/build/libspecbunch.dylib --output /external/build/grid-tests.json` checks automatic expansion in both directions, both hard caps and their edge flags, and nearly equal floating-point initial endpoints through the C ABI. The independent reference uses the same documented search bounds and spacing.

`cauchy_boundary.py --output /external/build/cauchy.json` separately checks zero Gaussian variance against the analytic Cauchy reference.

The standard `tests/xlang/generate.py` driver includes these four cases automatically.

`generate.py --library /external/build/libspecbunch.dylib --root /external/build/xlang/fixtures` creates four compact fixtures through an independent ctypes C-ABI adapter: exact observations, weighted observations with a nonzero kink, mixed Voigt friction, and a joint mass fit. The existing `tests/xlang/run_python.py`, `run_r.R`, `run_julia.jl`, `run_stata.py`, and `tests/no_stata/run.py` consume this manifest and compare every named scalar, histogram column, curve column and penalty-path column at 1/2/4/8/default threads. Store generated manifests, samples, environments and receipts under the external build root.

No plots are produced. Error tolerances apply to the independent numerical comparison, not to statistical uncertainty.
