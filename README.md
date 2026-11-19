# specbunch

Bunching mass and buncher mean estimators in Python, R, Julia and Stata, sharing a deterministic C++ core.

Methods: `cosine` (default), `joint`. Friction: `gaussian` (default), `voigt`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

See the repository `docs/REFERENCE.md` for every public option, field, flag and refusal reason, and `VENDOR.md` for source provenance. Build artifacts must be placed outside the source tree.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves. Each language wrapper draws default plots: `mass` shows the empirical transform with its fitted curve and selected window, and `npmle` shows the fitted distributions; pass `plot=False` (Stata `noplot`) to skip them.

Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.

Postestimation also provides `nonbuncher_mean` (a local characteristic mean with the buncher mean fixed) and ETI by the Saez trapezoid with continuous completion. See the public reference for language-specific syntax, conditional SE interpretation and flags.
