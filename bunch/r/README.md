# specbunch

Bunching mass and buncher mean estimators in Python, R, Julia and Stata, sharing a deterministic C++ core.

Methods: `cosine` (default), `joint`. Friction: `gaussian` (default), `voigt`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

See the repository `docs/REFERENCE.md` for every public option, field, flag and refusal reason, and `VENDOR.md` for source provenance. Build artifacts must be placed outside the source tree.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves.

## Plotting

By default `npmle` draws two views of the fit next to the data. The density view shows the data histogram (observed share divided by bin width), the fitted observed density after friction (bunchers included) and the fitted latent density with the bunching atom excluded, with a note giving the fixed B and friction and a thin dotted line at the kink. The CDF view shows the empirical CDF, the fitted observed CDF and the fitted latent CDF, which jumps by B at the kink. The horizontal range shows the 0.1% to 99.9% empirical quantiles; this affects the display only, never the fit.

`npmle(..., plot=TRUE)` is the default and draws with base graphics (no extra dependency) in interactive sessions only, so `Rscript` runs do not open devices. An explicit `plot="both"`, `"pdf"` or `"cdf"` also draws in batch, for example into `png()`; `plot=FALSE` skips plotting. `plot(result, kind="cdf")` redraws a fit and returns it invisibly. `mass(..., npmle=TRUE)` draws the default view when interactive.

Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.
