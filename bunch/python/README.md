# specbunch

Bunching mass and buncher mean estimators in Python, R, Julia and Stata, sharing a deterministic C++ core.

Methods: `cosine` (default), `joint`. Friction: `gaussian` (default), `voigt`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

See the repository `docs/REFERENCE.md` for every public option, field, flag and refusal reason, and `VENDOR.md` for source provenance. Build artifacts must be placed outside the source tree.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves.

## Plotting

By default `npmle` draws two views of the fit next to the data. The density view shows the data histogram (observed share divided by bin width), the fitted observed density after friction (bunchers included) and the fitted latent density with the bunching atom excluded, with a note giving the fixed B and friction and a thin dotted line at the kink. The CDF view shows the empirical CDF, the fitted observed CDF and the fitted latent CDF, which jumps by B at the kink. The horizontal range shows the 0.1% to 99.9% empirical quantiles; this affects the display only, never the fit.

`npmle(..., plot="both")` is the default; `plot="pdf"` or `"cdf"` draws one view and `plot=False` (or `None`) skips plotting. Plotting uses matplotlib, an optional dependency: `pip install "specbunch[plot]"`. Without matplotlib, `npmle` returns its results and issues one warning. The figure is stored in `result.figure` (a matplotlib `Figure`; save it with `result.figure.savefig("npmle.png")`), and `result.plot(kind="cdf")` draws a new one. The package never calls `plt.show()`. When pyplot is already imported or IPython is running, the figure is managed by pyplot, so it displays in notebooks or after `plt.show()`; otherwise it is a standalone `Figure`. Pass `plot=False` inside simulation loops.

`mass(..., plot=True)` is also the default: it draws the empirical cosine transform with a ±2 SE band, the fitted curve over the whole grid and the shaded selected window (joint fits add a t·S(t) panel) into `result.figure`, on scaled frequency t·σ_pilot. `plot="raw"` uses raw t and `plot=False` skips it; `result.plot(frequency="raw")` draws a new figure. With `npmle=True`, the NPMLE figure is `result.npmle.figure`.

Run local tests with `python -m pytest tests`; expensive checks marked `slow` are skipped by default. Use `--runslow` (or `SPECBUNCH_RUN_SLOW=1`) to include them. Full wheel CI enables them on pushes to `main` and manual dispatches.

Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.
