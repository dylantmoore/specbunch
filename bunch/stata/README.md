# specbunch

Bunching mass and buncher mean estimators in Python, R, Julia and Stata, sharing a deterministic C++ core.

Methods: `cosine` (default), `joint`. Friction: `gaussian` (default), `voigt`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

See the repository `docs/REFERENCE.md` for every public option, field, flag and refusal reason, and `VENDOR.md` for source provenance. Build artifacts must be placed outside the source tree.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves.

## Plotting

By default `npmle` draws two views of the fit next to the data. The density view shows the data histogram (observed share divided by bin width), the fitted observed density after friction (bunchers included) and the fitted latent density with the bunching atom excluded, with a note giving the fixed B and friction and a thin dotted line at the kink. The CDF view shows the empirical CDF, the fitted observed CDF and the fitted latent CDF, which jumps by B at the kink. The horizontal range shows the 0.1% to 99.9% empirical quantiles; this affects the display only, never the fit.

`specbunch npmle` draws a `twoway` graph by default, combining the two views with `graph combine`. Use `plot(pdf)`, `plot(cdf)` or `plot(both)` to choose the view and `noplot` to skip it; `name()` and `saving()` are passed to the final graph, for example `name(npmle, replace)` or `saving(npmle.gph, replace)`. Export with `graph export npmle.png, name(npmle)`. Drawing uses `preserve`/`restore`, so the data and `e(sample)` are unchanged.

`specbunch mass` also draws by default: the empirical cosine transform with a ±2 SE band, the fitted curve over the whole grid and the shaded selected window, on scaled frequency t·σ_pilot (`plot(scaled)`). `plot(raw)` uses raw t and `noplot` skips it; `name()` and `saving()` work as for `npmle`. Joint fits combine a cosine panel and a t·S(t) panel. Give the mass and NPMLE graphs different names to keep both.

Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.

Postestimation also provides `nonbuncher_mean` (a local characteristic mean with the buncher mean fixed) and ETI by the Saez trapezoid with continuous completion. See the public reference for language-specific syntax, conditional SE interpretation and flags.
