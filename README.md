# specbunch

Bunching mass and buncher mean estimators in Python, R, Julia and Stata, sharing a deterministic C++ core.

Methods: `gaussian`, `voigt`, `gaussian_jump`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

See the repository `docs/REFERENCE.md` for every public option, field, flag and refusal reason, and `VENDOR.md` for source provenance. Build artifacts must be placed outside the source tree.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves. It returns data for plotting and contains no plotting code.
