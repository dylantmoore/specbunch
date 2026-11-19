# specbunch

Bunching mass and buncher mean estimators in Python, R, Julia and Stata, sharing a deterministic C++ core.

## Installation

Nothing is on CRAN yet; this compiles the package from source (macOS: Xcode Command Line Tools; Windows: Rtools; Linux: `r-base-dev`):

```r
install.packages("remotes")
remotes::install_github("dylantmoore/specbunch", subdir = "bunch/r")
```

Pre-release research software: read [Status and testing](https://github.com/dylantmoore/specbunch#status-and-testing) before relying on it. Cite: Moore, Dylan (2026). "An Application of Euler's Formula to Bunching Methodology: Addressing Frictions via Spectral Deconvolution." Working paper.

Methods: `cosine` (default), `joint`. Friction: `gaussian` (default), `voigt`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

See the repository `docs/REFERENCE.md` for every public option, field, flag and refusal reason, and `VENDOR.md` for source provenance. Build artifacts must be placed outside the source tree.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves.

## Plotting

By default `npmle` draws two views of the fit next to the data. The density view shows the data histogram (observed share divided by bin width), the fitted observed density after friction (bunchers included) and the fitted latent density with the bunching atom excluded, with a note giving the fixed B and friction and a thin dotted line at the kink. The CDF view shows the empirical CDF, the fitted observed CDF and the fitted latent CDF, which jumps by B at the kink. The horizontal range shows the 0.1% to 99.9% empirical quantiles; this affects the display only, never the fit.

`npmle(..., plot=TRUE)` is the default and draws with base graphics (no extra dependency) in interactive sessions only, so `Rscript` runs do not open devices. An explicit `plot="both"`, `"pdf"` or `"cdf"` also draws in batch, for example into `png()`; `plot=FALSE` skips plotting. `plot(result, kind="cdf")` redraws a fit and returns it invisibly. `mass(..., npmle=TRUE)` draws the default view when interactive.

`mass(..., plot=TRUE)` follows the same rules for the transform fit: the empirical cosine transform with a ±2 SE band, the fitted curve over the whole grid and the shaded selected window (joint fits add a t·S(t) panel), on scaled frequency t·σ_pilot. `plot="raw"` uses raw t and `plot=FALSE` skips it; `plot(result, frequency="raw")` redraws. With `npmle=TRUE`, interactive sessions draw the transform fit and then the NPMLE plot.

Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.

## Licence and citing

Our code is MIT licensed, Copyright (c) 2026 Dylan Moore. Third-party components
retain their own licences; see THIRD_PARTY_NOTICES.md and the supplied licence texts.
Please cite Moore, Dylan (2026). *An Application of Euler's Formula to Bunching Methodology: Addressing Frictions via Spectral Deconvolution*. Working paper.
This paper title is a placeholder and will be updated. See CITATION.cff.
