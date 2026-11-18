# specbunch

Bunching mass and buncher mean estimators in Python, R, Julia and Stata, sharing a deterministic C++ core.

Methods: `gaussian`, `voigt`, `gaussian_jump`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

See the repository `docs/REFERENCE.md` for every public option, field, flag and refusal reason, and `VENDOR.md` for source provenance. Build artifacts must be placed outside the source tree.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves.

## Plotting

By default `npmle` draws two views of the fit next to the data. The density view shows the data histogram (observed share divided by bin width), the fitted observed density after friction (bunchers included) and the fitted latent density with the bunching atom excluded, with a note giving the fixed B and friction and a thin dotted line at the kink. The CDF view shows the empirical CDF, the fitted observed CDF and the fitted latent CDF, which jumps by B at the kink. The horizontal range shows the 0.1% to 99.9% empirical quantiles; this affects the display only, never the fit.

Plotting is a package extension for Plots.jl: load it with `using Plots`. `npmle(...; plot=:both)` is the default; `:pdf` or `:cdf` draws one view and `plot=false` skips plotting. With Plots loaded the figure is stored in `result.figure` and displayed in interactive sessions, and `plot(result; kind=:cdf)` redraws via the recipe (Plots attributes such as `xlabel` can be overridden). Without Plots the default request is a no-op with a one-time `@info`. `test/plots.jl` is a smoke test for an environment that has Plots.
