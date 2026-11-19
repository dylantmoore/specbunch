---
name: specbunch
description: Use when estimating bunching at a kink (tax kink, benefit kink, threshold in a schedule) where optimization frictions spread bunchers into a diffuse hump rather than an exact spike, using the specbunch package in Python, R, Julia or Stata. Covers the bunching mass B and friction scale, the mean characteristic of bunchers (buncher_mean) and of nonbunchers at the kink (nonbuncher_mean), the elasticity of taxable income (eti), and an NPMLE counterfactual (latent) density with the atom removed (npmle); also how to read refusals, tiers, flags and the windows table.
---

# specbunch

**What it estimates.** The share B of the population bunching at a kink, when each observation is measured with a symmetric friction (Gaussian, or Gaussian plus Cauchy "Voigt"); then the buncher and nonbuncher characteristic means, the ETI, and the latent density. **Key idea:** friction multiplies the characteristic function by the friction's characteristic function phi(t), so the atom appears as `B * phi(t)` in the empirical cosine transform at frequencies where the smooth density has decayed; B and the friction scale are fitted there by weighted nonlinear least squares over ranked frequency windows, with checks that the implied distribution is coherent.

**Install:** see the README (`Installation`). Everything builds from source today (C++20 compiler needed). Reference for every field: `docs/REFERENCE.md`.

## API

| Task | Python (`import specbunch`) | R (`library(specbunch)`) | Julia (`using SpecBunch`) | Stata |
|---|---|---|---|---|
| Bunching mass | `specbunch.mass(z, kink=0, weights=None, method="cosine", friction="gaussian", policy="tiered", seed=1, threads=None, check_all_windows=False, npmle=False, plot=True)` | `mass(z, kink=0, weights=NULL, method="cosine", friction="gaussian", policy="tiered", seed=1, threads=NULL, check_all_windows=FALSE, npmle=FALSE, plot=TRUE)` | `mass(z; kink=0.0, weights=nothing, method=:cosine, friction=:gaussian, policy=:tiered, seed=1, threads=nothing, check_all_windows=false, npmle=false, plot=true)` | `specbunch mass z [pw=w], kink(0) method(cosine) friction(gaussian) policy(tiered) seed(1) [checkallwindows noplot plot(raw)]` |
| Buncher mean of x | `specbunch.buncher_mean(z, x, kink=0, weights=None, seed=1, draws=499)` | `buncher_mean(z, x, kink=0, weights=NULL, seed=1, draws=499)` | `buncher_mean(z, x; kink=0.0, weights=nothing, seed=1, draws=499)` | `specbunch buncher z x [pw=w], kink(0) draws(499)` |
| Nonbuncher mean at kink | `specbunch.nonbuncher_mean(z, x, mass_fit, buncher_mean_fit, weights=None, bandwidth=.05, degree=1)` (or `mu_b=`, `se_mu_b=` instead of a buncher fit) | `nonbuncher_mean(z, x, mass_fit, buncher_mean_fit, weights=NULL, bandwidth=.05, degree=1)` | `nonbuncher_mean(z, x, mass_fit, buncher_mean_fit; weights=nothing, bandwidth=.05, degree=1)` | `specbunch nonbuncher_mean z x [pw=w], massfit(m) buncherfit(b) bandwidth(.05) degree(1)` (or `mub(#) semub(#)`) |
| ETI | `specbunch.eti.estimate(fit, z, tax_below=t0, tax_above=t1, weights=None)` (or `net_below=`, `net_above=`) | `eti(fit, z, tax_below=t0, tax_above=t1, weights=NULL)` | `SpecBunch.ETI.estimate(fit, z; tax_below=t0, tax_above=t1, weights=nothing)` | `specbunch eti z [pw=w], massfit(m) taxbelow(#) taxabove(#)` (or `netbelow(#) netabove(#)`) |
| NPMLE latent density | `specbunch.npmle(fit, z, weights=None, bins=None, plot="both")`, or `mass(..., npmle=True)` | `npmle(fit, z, weights=NULL, bins=NULL, plot=TRUE)` | `npmle(fit, z; weights=nothing, bins=nothing, plot=:both)` | `specbunch npmle [, plot(both) bins(matname) noplot]` after `specbunch mass` |
| Candidate windows | `fit.windows` (pandas DataFrame; dict of arrays without pandas) | `fit$windows` | `fit.windows` (NamedTuple of vectors) | `specbunch windows` (`e(windows)`) |

Options: `method` is `cosine` (default; cosine transform only) or `joint` (also fits the sine transform, allowing a density jump at the kink; returns `jump`). `friction` is `gaussian` (returns `sigma`) or `voigt` (returns Cauchy scale `a`, Gaussian variance `s2`, and `sigma_eq`). `policy` is `tiered` (default) or `first`. In Stata, store first-stage fits with `estimates store` and pass their names to `massfit()` / `buncherfit()`. Plots: Python needs the `[plot]` extra; Julia needs `using Plots`; pass `plot=False` / `noplot` in loops.

## Reading results

- **`status` / `refused` / `reason`.** `status` 0 is an answer. Refusals: 1 `insufficient_data` (too few usable rows), 2 `no_frequency_band` (no pilot scale or visible band), 3 `no_coherent_window` (no window gives 0 < B < 1 with positive friction). A refusal is a result, not a bug: report it rather than trimming or reshaping the sample until it goes away.
- **`tier`** (tiered policy): 1 all checks pass; 2 fit test passes but a distribution-shape check flags; 3 shape checks pass but `poor_fit`; 4 any coherent window. Report the tier.
- **`flags`** (a coherent estimate is returned even when flagged):
  - `poor_fit`: the window's fit-test p-value (`fit_p`) is below .10.
  - `cdf_inconsistent`: the implied distribution violates the DKW band around the empirical CDF.
  - `moments_inconsistent`: violates simultaneous harmonic-moment bounds.
  - `no_valid_distribution`: an LP certificate shows no distribution has the fitted moments.
  - `distribution_check_inconclusive`: the LP budget ran out; this is not a pass and not a failure.
  - `residual_cdf_inconsistent`: the probit residual-CDF screen rejects.
- **`windows` table:** one row per candidate frequency window in rank order, with B, sigma, SEs, `fit_p`, coherence and each check as `pass`, `fail` or `not_evaluated` (never read `not_evaluated` as pass; `check_all_windows=True` evaluates all without changing the answer), and `selected`. Use it as a sensitivity table: does B move much across coherent, passing windows?
- **Default plots** are diagnostics: the transform plot shows the empirical cosine transform with ±2 SE, the fitted `B * phi(t)` curve and the selected window shaded (joint adds a t·S(t) panel); NPMLE plots compare data, fitted observed and latent density/CDF.
- **Standard errors** (`se_log_B`, `se_log_sigma`, buncher `se`, nonbuncher `se`) are model-based, conditional on the selected window, and **not confidence intervals**. ETI has no SE.
- Buncher mean: `mu` (with `pointwise_median`, `se`, `poor_fit`), bands in its `windows` table. Nonbuncher mean: `m_NB`, `pi_kink`, flags `outside_outcome_range`, `density_floor_used`, `sparse_bins`, `quadrature_inconclusive`. ETI: `saez.e` (headline) and `continuous.e` with `.d`, each with a `status`; flags `negative_shoulder`, `zero_shoulder`.

## Recommended workflow

1. Build `z` in logs centred at the kink and keep the full sample.
2. Fit `mass` with `friction="gaussian"` and `friction="voigt"`; report both as a sensitivity pair (nonbuncher mean and ETI depend on the friction law).
3. Check `status`, `tier`, `flags`, look at the transform plot and the `windows` table before interpreting B.
4. Then `buncher_mean`, `nonbuncher_mean`, `eti` or `npmle` from a coherent mass fit, with the same sample, weights and kink.
5. Report the platform (OS, CPU), package version (`version()`) and seed; results are bitwise reproducible on one platform at any thread count.

## Pitfalls

- **Units.** `z` is log income minus log kink (or log income with `kink=log(K)`). The friction scale is then in log points, and `eti` requires log income in the same coordinate as the mass fit. The package never logs `z` for you.
- **Full sample, no trimming or windowing.** The method uses the whole distribution; do not cut to a band around the kink or drop tails. Missing, nonfinite and nonpositive-weight rows are dropped automatically.
- **Weights** are positive observation (sampling) weights; only relative size matters (rescaling all weights changes nothing beyond rounding); the effective sample size is Kish's `(sum w)^2 / sum w^2`. Stata accepts `pweight` and `fweight` but treats both this way, so an `fweight` is not expanded into repeated observations. Pass the same weights to every stage.
- **`method="joint"` is unweighted only**; weights raise an error. With Voigt friction only the `poor_fit` screen runs (shape checks stay `not_evaluated`).
- `draws` (buncher mean) must be at least 19. `nonbuncher_mean` `bandwidth` is a declared sensitivity choice, not optimized; `degree=2` is a sensitivity check.
- Estimators and defaults may change before the paper is final; pin the commit you used.

## Citation

Moore, Dylan (2026). "An Application of Euler's Formula to Bunching Methodology: Addressing Frictions via Spectral Deconvolution." Working paper.

The paper title is a placeholder and will be updated; see CITATION.cff.
