# specbunch

Bunching mass and buncher mean estimators in Python, R, Julia and Stata, sharing a deterministic C++ core.

## Installation

The binary package `SpecBunch_jll` is not registered yet, so `SpecBunch` is not in the General registry. For now, build the core library with CMake and develop the package against a local stand-in; the copy-paste steps are in the [main README](https://github.com/dylantmoore/specbunch#julia-19). Pre-release research software: read [Status and testing](https://github.com/dylantmoore/specbunch#status-and-testing) before relying on it. Cite: Moore, Dylan (2026). "Spectral Deconvolution: Applications to Bunching and Regression Discontinuity Designs." Working paper.

Methods: `cosine` (default), `joint`. Friction: `gaussian` (default), `voigt`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

See the repository `docs/REFERENCE.md` for every public option, field, flag and refusal reason, and `VENDOR.md` for source provenance. Build artifacts must be placed outside the source tree.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves.

## Plotting

By default `npmle` draws two views of the fit next to the data. The density view shows the data histogram (observed share divided by bin width), the fitted observed density after friction (bunchers included) and the fitted latent density with the bunching atom excluded, with a note giving the fixed B and friction and a thin dotted line at the kink. The CDF view shows the empirical CDF, the fitted observed CDF and the fitted latent CDF, which jumps by B at the kink. The horizontal range shows the 0.1% to 99.9% empirical quantiles; this affects the display only, never the fit.

Plotting is a package extension for Plots.jl: load it with `using Plots`. `npmle(...; plot=:both)` is the default; `:pdf` or `:cdf` draws one view and `plot=false` skips plotting. With Plots loaded the figure is stored in `result.figure` and displayed in interactive sessions, and `plot(result; kind=:cdf)` redraws via the recipe (Plots attributes such as `xlabel` can be overridden). Without Plots the default request is a no-op with a one-time `@info`. `mass(...; plot=true)` is also the default: with Plots loaded it stores the transform fit in `result.figure` (empirical cosine transform with a ±2 SE band, the fitted curve over the whole grid and the shaded selected window; joint fits add a t·S(t) panel), on scaled frequency t·σ_pilot. `plot=:raw` uses raw t and `plot=false` skips it; `plot(result; frequency=:raw)` redraws. With `npmle=true`, the NPMLE figure is `result.npmle.figure`. `test/plots.jl` is a smoke test for an environment that has Plots.

Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.

Postestimation also provides `nonbuncher_mean` (a local characteristic mean with the buncher mean fixed) and ETI by the Saez trapezoid with continuous completion. See the public reference for language-specific syntax, conditional SE interpretation and flags.
