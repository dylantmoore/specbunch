# Public reference

`specbunch` estimates a bunching mass and a buncher outcome mean using frequency-domain fits. Python, R, Julia and Stata call the same deterministic C++ core. Mass estimation draws a default transform-fit plot, and NPMLE postestimation provides default density and CDF plots, through the language wrappers.

## Mass

The headline estimator is weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits the sine transform. The companions are log-linear OLS and LAD on the reference band.

Python: `mass(z, kink=0, weights=None, method="cosine", friction="gaussian", policy="tiered", threads=None, seed=1, check_all_windows=False, lp_iterations=5000, lp_retry_iterations=31250, npmle=False, plot=True)`.
R uses the same names with `NULL` and `FALSE`; Julia uses keyword arguments and method/friction/policy symbols. Stata uses `specbunch mass z [if] [in] [pw=weight], method(cosine) friction(gaussian) policy(tiered) [plot(scaled|raw) noplot name() saving()]`.

| Option | Meaning |
|---|---|
| z | One-dimensional observed running variable. |
| kink | Location of the bunching atom, subtracted before fitting; default zero. |
| weights | Optional positive observation weights; rows with nonfinite inputs or nonpositive weights are removed together. |
| method | `cosine` (default) fits only the empirical cosine transform; `joint` jointly fits cosine and sine, allowing a density jump at the kink. |
| friction | `gaussian` (default) uses Gaussian friction; `voigt` uses symmetric Gaussian plus Cauchy friction. |
| policy | `tiered` selects by the tiers below; `first` selects the highest-ranked coherent window and reports failed checks as flags. |
| threads | Core thread count; omitted, zero or negative uses `OMP_NUM_THREADS` when set, otherwise the core default. |
| seed | Nonnegative integer bootstrap seed, with neutral default 1. |
| check_all_windows | Evaluate distribution-shape checks on every coherent window when true; otherwise evaluate only as needed during selection. |
| lp_iterations | Initial LP simplex iteration budget, default 5000. |
| lp_retry_iterations | Retry LP simplex iteration budget, default 31250. |
| npmle | Attach spline postestimation with its default options when true and the mass estimate is coherent; default false. |
| plot | Draw the transform fit by default, on scaled frequency; `raw` uses raw frequency and false skips drawing. See [Transform fit plot](#transform-fit-plot). |

Stata spells the window-check and LP options `checkallwindows`, `lpiterations()` and `lpretryiterations()`; `kink()`, `threads()` and `seed()` have the meanings above.

`method="joint"` requires unweighted observations. Invalid argument types or unsupported options raise errors; inadequate statistical information returns a refused result. `loglinear_ols` and `loglinear_lad` report log-linear OLS and LAD estimates on the fixed reference frequency band alongside each mass result.

The default rule selects the highest-ranked coherent window in the first nonempty tier: (1) all checks pass; (2) fit test passes, with distribution-shape flags allowed; (3) distribution-shape checks pass, with a poor-fit flag allowed; (4) any coherent window. An inconclusive distribution check is distinguished from a demonstrated shape failure. The selected tier is returned as `tier` (1–4). The `first` rule does not use checks to select the window.

## Coherence, flags and refusals

| Field | Meaning |
|---|---|
| fit_converged | The numerical window fit converged. |
| coherent_values | The fit converged with 0 < B < 1 and positive friction scale. |
| cdf_inconsistent | The implied distribution violates the DKW CDF band. |
| moments_inconsistent | The implied distribution violates simultaneous harmonic-moment bounds. |
| no_valid_distribution | An LP certificate excludes every distribution with the specified moments. |
| distribution_check_inconclusive | The LP budget ended without a witness or exclusion certificate. |
| residual_cdf_inconsistent | The probit residual-CDF screen rejects. |
| poor_fit | The selected window's fit-test p-value is below .10. |
| insufficient_data | Too few usable observations remain to fit the estimator. |
| no_frequency_band | The pilot scale or frequency visibility calculation produces no usable band. |
| no_coherent_window | No candidate window or band yields a coherent estimate. |

The first two fields define coherence and are not warning flags. A coherent answer is returned even when checks fail. `status` is 0 for an answer, 1 for insufficient data, 2 for no frequency band and 3 for no coherent window; `refused` is the corresponding Boolean and `reason` is empty for an answer. Scalar check fields are 1 when flagged and 0 otherwise; `flags` is a list of names in language wrappers and the raw integer bit mask uses the six flags above in order, starting at bit zero.

## Returned mass fields and windows

`B` is the estimated atom share; `sigma` is Gaussian friction standard deviation; `se_log_B` and `se_log_sigma` are model-based standard errors on log scales; `fit_p` is the window fit-test p-value. Voigt additionally returns Cauchy scale `a`, Gaussian variance `s2`, equivalent scale `sigma_eq` and its `se_log_sigma_eq`. For Voigt, `se_log_sigma` is the retained alias of `se_log_sigma_eq`, not a standard error for the Gaussian component alone; its coherence scale is `sigma_eq`, so pure Cauchy friction is allowed. The joint fit additionally returns `jump`, the coefficient of the sine-transform term phi(t)/t associated with a right-minus-left density jump. Standard errors are conditional on the selected fit and are not confidence intervals.

`rank` identifies the selected candidate in rank order; `window_lo` and `window_hi` are its frequency endpoints and the `_scaled` versions multiply these by `sigma_pilot`. `sigma_pilot` is the pilot friction scale, `n` the retained row count, `n_input` the input row count, `n_eff` the effective sample size when provided, `W` the total weight, `weighted` the weight indicator, `seed` the seed used and `threads_used` the actual thread count. `method`, `friction`, `estimator`, `policy`, `se_note`, `flags_str` and `check_outcomes` describe the estimate and raw check encoding. Simple results have `B`, `sigma`, `status`, `m` (frequency count) and `ndrop` (discarded frequency count).

The default `windows` table has one row per candidate in rank order, with `rank`, `window_lo`, `window_hi`, `window_lo_scaled`, `window_hi_scaled`, `n_points`, `B`, `sigma`, `se_log_B`, `se_log_sigma`, `fit_p`, `fit_converged`, `coherent_values`, `coherent`, the six check columns above, and `selected`. Voigt adds `a`, `s2`, `sigma_eq`; the joint method adds `jump`. `coherent` and `selected` are 0/1 indicators. Each check cell is `pass`, `fail` or `not_evaluated`; raw C and Stata encode these as 1, 0 and -1. Unevaluated checks must not be interpreted as passing. Enabling all-window checks does not change the selected answer.

Python returns a pandas DataFrame when pandas is installed, otherwise a dictionary of NumPy arrays; R returns a data.frame; Julia returns a NamedTuple of vectors. Stata stores `e(windows)` with column names; `specbunch windows` lists it. `raw` preserves the named C outputs in Python, R and Julia, including numeric `windows.*` columns.

The simple estimators share one reference band, chosen before either log-linear fit using the pilot visibility cutoff and a levels-fit screen; they discard nonpositive cosine coefficients and do not select separate bands. Their `status` is 0 when available, 1 when the pilot scale is unavailable, 2 when the pilot mass is nonpositive, 3 when there is no visibility cutoff, 4 when there are no reference-band candidates, 5 when no valid reference fit exists, 6 when the log-linear fit is invalid, and 7 when fewer than three positive cosine coefficients remain. Their `m` counts reference-band frequencies before discarding coefficients, and `ndrop` counts discarded coefficients. These side-estimate status codes do not control refusal of the main mass estimate.

## Buncher mean

Python: `buncher_mean(z, x, kink=0, weights=None, threads=None, seed=1, draws=499)`; R and Julia use the same option names. Stata: `specbunch buncher z x, draws(499)`. Inputs, weights, kink, seed and threads have the meanings above; `x` is the outcome and `draws` is the fit-test bootstrap count (at least 19).

The estimator selects a ranked band with positive denominator and a mean within the retained outcome range. A band that fails only its fit test remains usable with `poor_fit`; if no coherent band passes the fit test, the first coherent band is selected. `mu` is the cosine-ratio mean, `pointwise_median` is the median frequency-specific ratio, and `se` is its model-based standard error. `band_lo` and `band_hi` are frequency endpoints.

The default band table contains `rank`, `band_lo`, `band_hi`, `band_lo_scaled`, `band_hi_scaled`, `mu`, `pointwise_median`, `se`, `fit_p`, `coherent`, `poor_fit`, `flags` and `selected`. Scaled endpoints multiply by the pilot scale. The numeric Stata table uses `poor_fit` as the flag indicator. `band_rule`, `p0_rule` and complex-ratio estimates are not part of the public API.

## Reproducibility and provenance

Pass an explicit seed to reproduce a particular fit-test realization. The same build produces bitwise-identical outputs across language wrappers and thread counts. Dependency and source provenance is recorded in VENDOR.md.

## Spline postestimation

Python: `npmle(fit, z, weights=None, bins=None, *, plot="both", threads=None, core_intervals=12, grid_points=1001, exact_threshold=1000, max_iterations=200, lambda_min=.01, lambda_max=10000, lambda_points=15, refinement_points=9, histogram_width=None, quadrature_tolerance=1e-9, optimality_tolerance=1e-7)`.
R and Julia expose the same names with their usual null and keyword conventions. `mass(..., npmle=True)` in Python, `npmle=TRUE` in R, and `npmle=true` in Julia attach the default postestimation result; a refused mass result remains a refusal without a spline fit. Stata uses `specbunch npmle` after `specbunch mass`.

The latent model is an atom of mass B at the kink plus a smooth density of total mass 1−B. Each side has its own clamped cubic B-spline basis, allowing distinct density limits at the kink. B and the friction parameters are fixed at the supplied mass estimate. The `friction` argument determines the convolution for either method: `gaussian` uses a Gaussian; `voigt` convolves a Gaussian of variance `s2` with a Cauchy of scale `a`, including the pure-Cauchy boundary.

The negative log-likelihood is convex in the nonnegative spline coefficients subject to their integral being 1−B. The objective adds `(lambda/2) theta' D'D theta` independently on each side, where D takes second differences of the spline heights. A deterministic active-set Newton method solves this constrained problem to the reported KKT tolerance. There is no sample truncation: every retained observation contributes, including remote tails. Observation weights are rescaled to sum to the Kish effective sample size `(sum w)^2/sum(w^2)`; provide the same weights used for the mass estimate when reproducing its weighted sample.

For samples larger than `exact_threshold`, the likelihood uses occupied bins anchored at the kink, of default width `max(sigma,a)/20`, with exact model bin probabilities rather than density values at bin centers. Empty bins need not be stored because they have zero multinomial counts. Smaller samples use the exact observed densities. `loglik` is the unpenalized sum of effective counts times log bin probabilities, or log observed densities for exact observations; it excludes multinomial combinatorial constants. These two likelihood conventions are not directly comparable across different histogram widths.

BIC is `−2 loglik + edf log(n_eff)`, with `edf = tr((F'WF + P)^{-1} F'WF)`, using the observed likelihood information and the second-difference penalty. The trace includes all spline coefficients and follows this stated convention even when some coefficients lie on the nonnegativity boundary. A QR factorization of the augmented weighted design evaluates the trace without forming an inverse; `edf_ridge` is zero. The initial lambda grid contains 15 logarithmically spaced values from .01 to 10000. If its BIC minimum is at an endpoint, the search extends in that direction using the initial log spacing (at least a quarter decade per extension) until the minimum is interior or reaches a cap of 1e-6 or 1e10. An interior minimum is refined with nine log-spaced values between its coarse neighbors; the grid-edge flag identifies a selected capped endpoint. Explicit initial endpoints beyond those caps remain endpoints.

### Main options

| Option | Meaning |
|---|---|
| fit | A coherent `mass` result supplying B, friction and kink. |
| z | The complete observed running-variable sample for postestimation. |
| weights | Optional positive observation weights, filtered jointly with z. |
| bins | Optional finite, strictly increasing histogram output edges in the original z units; these never restrict the likelihood sample. |
| plot | Draw density and CDF views by default; language-specific controls are described below. |

### Advanced (solver settings)

| Option | Meaning |
|---|---|
| threads | Shared core thread count; reductions have a fixed order. |
| core_intervals | Number of equal knot intervals on each side within the central support, default 12; geometrically widening tail intervals cover every observation. |
| grid_points | Number of equally spaced curve points in the central range, default 1001; additional tail points and kink limits cover the full sample range. |
| exact_threshold | Largest retained sample size evaluated observation by observation, default 1000. |
| max_iterations | Newton iteration budget per lambda, default 200; an unsuccessful warm start is retried from uniform mixture weights with twice this budget. |
| lambda_min, lambda_max | Positive endpoints of the initial penalty grid, default .01 and 10000; the search extends in the initial log spacing (at least a quarter decade) until the BIC minimum is interior or reaches 1e-6 or 1e10, then refines an interior minimum (explicit initial endpoints outside those caps are retained). |
| lambda_points | Number of coarse logarithmically spaced penalty values, default 15. |
| refinement_points | Number of penalty values used to refine an interior minimum, default 9. |
| histogram_width | Positive likelihood-bin width in z units; null uses max(sigma,a)/20. |
| quadrature_tolerance | Target absolute error in a CDF or in max(sigma,a) times a PDF for each Voigt kernel evaluation, default 1e-9. |
| optimality_tolerance | Target normalized simplex KKT residual during iteration, default 1e-7; final acceptance allows ten times this target (1e-6 by default) to tolerate line-search roundoff without changing the iteration path. |

The central spline support reaches the .5% and 99.5% sample quantiles, with at least four friction scales on each side. Additional knots double their spacing until the support extends beyond every observation by at least four friction scales. This uses quantiles to place knots, not to discard data. Output histogram bins have central width approximately max(sigma,a)/4, with at most 4000 central intervals and widening tail bins; user bins replace only this output grid. Curves include both the point immediately below the kink and the kink itself when it lies within the sample range, so the latent CDF's atom jump is represented explicitly.

Gaussian spline convolution uses closed-form truncated normal moments. The Gaussian omission probability is bounded by `erfc(9/sqrt(2))`, about 2.3e-19, and is reported as `gaussian_tail_bound`; `gaussian_pdf_error_bound` reports the corresponding conservative density bound. Voigt convolution integrates the analytic Gaussian convolution against the Cauchy kernel by deterministic adaptive Gauss–Legendre quadrature. `quadrature_error` is the accumulated numerical error estimate, not a certified analytic bound on quadrature roundoff; output histogram probabilities use differences of CDF values, while Voigt likelihood probabilities integrate the positive Cauchy bin probability directly to avoid tail cancellation. The independent simulated-data reference uses analytic Cauchy spline integrals and high-order Gaussian expectation. Numerical integration accuracy is validated against that separate calculation.

| Scalar or field | Meaning |
|---|---|
| lambda | Selected roughness penalty; Python also exposes `lambda_` because lambda is a language keyword. |
| edf | Effective degrees of freedom used in BIC. |
| loglik | Unpenalized full-sample likelihood under the stated exact or histogram convention. |
| bic | Penalized model-selection score at the selected lambda. |
| lambda_at_grid_edge | One if the selected penalty is at a capped search endpoint; an initial endpoint alone does not trigger this flag. |
| latent_left, latent_right | Left and right limits of the smooth latent density, including the factor 1−B. |
| latent_jump | Right minus left latent density limit; this jump is weakly identified. |
| status | `converged`, `solver_not_converged`, or `quadrature_not_converged`. |
| status_code | Corresponding native integer code 0, 1 or 2. |
| B, sigma, a, kink | Fixed mass, Gaussian standard deviation, Cauchy scale and kink used for postestimation. |
| n, n_input, n_eff | Retained, input and effective sample sizes. |
| histogram_width | Width used for the likelihood histogram. |
| exact_observations | One when the likelihood used exact observation densities. |
| likelihood_rows | Number of exact observations or occupied likelihood bins. |
| basis_count | Total number of spline basis functions across both sides. |
| coefficients | Nonnegative spline heights, ordered from kink to tail on the left, then from kink to tail on the right. |
| knots_left, knots_right | Nonnegative knot distances from the kink for each side. |
| kkt | Normalized simplex optimality residual at the selected fit. |
| iterations | Iterations used by the selected Newton solve. |
| edf_ridge | Zero: the QR-based EDF calculation requires no diagonal stabilization. |
| unconverged_lambdas | Number of evaluated penalty values whose numerical solve did not meet the requested tolerance. |
| quadrature_error | Maximum accumulated scaled quadrature-error estimate across evaluated rows and output points. |
| gaussian_tail_bound, gaussian_pdf_error_bound | Analytic Gaussian truncation probability and density bounds, respectively. |
| threads_used | Actual core thread count. |
| likelihood | `exact_observations` or `full_sample_histogram`. |
| jump_note | Reminder that the latent density jump is weakly identified. |
| lambda_path | Table of all evaluated lambda values, loglik, edf, bic, converged indicator and kkt. |

The `histogram` table contains `bin_lo`, `bin_hi` (original-unit edges), `observed_share` (weighted empirical probability), `fitted_share` (fitted observed probability), `nonbuncher_share` ((1−B) times the smooth component after friction), `buncher_share` (B times friction probability), and `latent_counterfactual_density` (smooth latent density at the bin midpoint). Shares use the full-sample denominator, including when custom output bins cover only part of the sample. Fitted probabilities outside the output range are not renormalized away.

The `curves` table contains `grid` (original-unit positions), `latent_pdf` ((1−B) times the smooth latent density, excluding the atom), `latent_cdf` (right-continuous latent CDF, including the atom), `observed_pdf` and `observed_cdf` (fitted distribution after friction, including bunchers), and `empirical_cdf` (weighted empirical CDF at the same points). Python, R and Julia also expose these six columns directly on the NPMLE result. Python returns DataFrames when pandas is present and NumPy-array dictionaries otherwise; R returns data.frames; Julia returns NamedTuples of vectors.

Stata preserves the mass result and stores `e(npmle_histogram)`, `e(npmle_curves)` and `e(npmle_lambda_path)` with column names; scalars have the `e(npmle_)` prefix and `e(npmle_status)` is the status string. Its option spellings are `bins(matrix_name)`, `threads()`, `coreintervals()`, `gridpoints()`, `exactthreshold()`, `maxiterations()`, `lambdamin()`, `lambdamax()`, `lambdapoints()`, `refinepoints()`, `histogramwidth()`, `quadraturetolerance()` and `optimalitytolerance()`, with the same defaults. The bins matrix must be a row or column vector. Postestimation uses the retained estimation sample and original weight expression, so the data and referenced variables must still be present.

### Plotting

Python defaults to `plot="both"`; `"pdf"` or `"cdf"` selects one view and `False` (or `None`) skips drawing. Install the optional `specbunch[plot]` extra for matplotlib. The figure is stored in `result.figure`; `result.plot(kind="both", xlabel=None)` returns a new figure. Missing matplotlib produces a warning for automatic plotting and leaves `figure=None`; the estimator still returns its tables. The wrapper does not call `plt.show()`.

R defaults to `plot=TRUE` and draws with base graphics only in interactive sessions when the argument is omitted. Explicit `TRUE`, `"both"`, `"pdf"` or `"cdf"` also draws in batch; `FALSE` or `"none"` suppresses drawing. The S3 method `plot(result, kind="both")` redraws the result.

Julia defaults to `plot=:both`; `:pdf` and `:cdf` select one view and `false` or `nothing` suppresses drawing. Load `Plots` to activate the package extension: the figure is stored in `result.figure` and displayed in interactive sessions. Without Plots, automatic plotting is a no-op with a one-time informational message. `Plots.plot(result; kind=:cdf)` redraws the result.

Stata draws by default with `specbunch npmle, plot(both)`; `plot(pdf)` or `plot(cdf)` selects one view and `noplot` suppresses it. `name(name[, replace])` and `saving(filename[, replace])` apply to the final graph. Weights are inherited from the mass estimation sample, rather than supplied again.

Density plots compare histogram density, fitted observed density and latent smooth density excluding the atom; CDF plots compare empirical, fitted observed and latent CDFs, including the latent atom jump. The display uses the central 0.1%–99.9% empirical range, expanded to include the kink; this does not restrict estimation or the returned full-range curves.

## Reproducibility contract

On the same platform and build, results are bitwise identical at any thread count.
Across operating systems, the selected answer agrees on the test fixtures within
relative tolerance 1e-8 (absolute tolerance 1e-12 near zero). Status, refusal reason,
tier, flags, and selected window endpoints must agree exactly; the numerical
comparison covers selected B, sigma, mu (buncher mean), and tau (regression discontinuity),
where applicable. The tolerance is max(1e-12, 1e-8 * max(abs(actual), abs(reference))).

System math libraries differ in the last bit. Diagnostics such as fit p-values and
window ranks can differ slightly, and a sample on a knife edge between two windows
can select a different window. Report the platform with published results.

CI keeps same-platform thread comparisons strictly bitwise. Cross-platform wheel
and Stata-harness jobs fail only for the selected-answer differences above. All
other differences, including candidate-table lengths and unselected windows, are
reported in the job summary and a compressed JSON-lines artifact with every field,
fixture, actual/reference value, and relative difference. Selected-answer fixture
flips are counted explicitly. NPMLE curves have a looser reporting tolerance of
relative 1e-6 or absolute 1e-8 (the larger allowance); exceeding it is reported,
not a cross-platform failure. These tolerances do not change estimator numerics.

The C ABI retains version 1, existing option layouts and existing entry points. The additive `sb_mass(data, opts, method, friction, out)` entry point accepts `SB_METHOD_COSINE` or `SB_METHOD_JOINT` and `SB_FRICTION_GAUSSIAN` or `SB_FRICTION_VOIGT`.

## Transform fit data

Every mass result includes a numeric `transform` table (Python pandas DataFrame when available, otherwise a dict of NumPy arrays; R data.frame; Julia NamedTuple of vectors; Stata e(transform)). The default mass plot draws this table; see below.
The table uses 401 equally spaced frequencies, from 0.1/pilot below the lowest candidate endpoint (clipped at zero) to 0.2/pilot beyond the highest candidate endpoint. It is empty when no candidate frequency range exists. Refused results with a range retain empirical curves but have missing fitted curves.

| Column | Meaning |
|---|---|
| t | Frequency in the units of the centered running variable. |
| t_scaled | Frequency multiplied by sigma_pilot. |
| cosine | Empirical cosine transform at t. |
| cosine_se | Pointwise standard error of the empirical cosine transform from its own variance. |
| cosine_fit | B times the fitted friction transform, evaluated over the whole grid. |
| sine_scaled | For joint fits, t times the empirical sine transform. |
| sine_scaled_se | For joint fits, t times the pointwise standard error of the empirical sine transform. |
| sine_fit | For joint fits, jump times the fitted friction transform (the fitted t-scaled sine curve). |

`window_lo`, `window_hi` and their scaled versions identify the selected window; `n_points` gives its number of fitting frequencies. These frequencies are distinct from the display grid. Weighted pointwise variances use sum(w_i^2 (f_i - weighted_mean(f))^2) / sum(w_i)^2, with no degrees-of-freedom correction. Curves are diagnostic evaluations after selection and never feed back into the fit.

### Transform fit plot

By default `mass` draws the transform fit. The first panel shows the empirical cosine transform at each grid frequency with a pointwise ±2 SE band (`cosine` ± 2 `cosine_se`), the fitted curve `cosine_fit` over the whole grid and the selected window shaded between dashed lines at its endpoints. A panel title gives method, friction and the estimates (B and sigma; B, a and s2 for Voigt). Joint fits add a second panel with `sine_scaled` ± 2 `sine_scaled_se`, `sine_fit` and the same window, titled with the jump estimate. The horizontal axis is frequency times `sigma_pilot` by default, or raw frequency t on request. A refused result with a frequency range shows the empirical curves only; a result with an empty table draws nothing.

Python defaults to `plot=True` (the same as `"scaled"`); `"raw"` uses raw frequency and `False`, `None` or `"none"` skips drawing. The figure is stored in `result.figure` and `result.plot(frequency="scaled")` returns a new one. The optional matplotlib extra, missing-matplotlib warning and no-`plt.show()` rules match NPMLE plotting.

R defaults to `plot=TRUE` and draws with base graphics only in interactive sessions when the argument is omitted. Explicit `TRUE`, `"scaled"` or `"raw"` also draws in batch; `FALSE` or `"none"` suppresses drawing. The S3 method `plot(result, frequency="scaled")` redraws the result and returns it invisibly.

Julia defaults to `plot=true` (the same as `:scaled`); `:raw` uses raw frequency and `false` or `nothing` suppresses drawing. With Plots loaded the figure is stored in `result.figure` and displayed in interactive sessions; without Plots the request is a no-op with a one-time informational message. `Plots.plot(result; frequency=:raw)` redraws the result.

Stata draws by default with `specbunch mass z, plot(scaled)`; `plot(raw)` uses raw frequency and `noplot` suppresses drawing. `name(name[, replace])` and `saving(filename[, replace])` apply to the final graph; joint fits combine the two panels with `graph combine`. Drawing uses `preserve`/`restore`.

With `npmle=True` (Python, R, Julia), the NPMLE figure is separate: `result.figure` holds the transform fit and `result.npmle.figure` holds the NPMLE plot. R draws the two plots in turn on the current device in interactive sessions; `plot(result)` and `plot(result$npmle)` redraw either. In Stata, `specbunch npmle` draws its own graph, so give the two graphs different `name()` values to keep both.

## Joint transform fit with Voigt friction

All four method/friction combinations are supported. For symmetric Voigt friction, `phi(t) = exp(-a*abs(t) - s2*t*t/2)`; the joint model fits `B*phi(t)` to the empirical cosine transform and `jump*phi(t)/t` to the empirical sine transform. Convolution by symmetric friction multiplies both Fourier components by the same real characteristic function, so the jump term carries the same phi. The covariance weights are fixed from the empirical cosine/sine covariance at each frequency.

The joint Voigt solver profiles the two linear coefficients (B and the scaled jump) exactly, then uses constrained Gauss–Newton steps with line search for nonnegative a and s2, checking both friction boundaries. Ranked windows use the joint Gaussian frequency-grid and visibility-extension design: at least 11 frequencies, a pilot-scaled spacing of 0.05 from 0.2 to an initial upper bound of 4, extended up to 8 when signal remains visible. Its fit test uses 499 deterministic sign-multiplier draws, refitting the joint model and comparing the maximum standardized cosine/sine residual. Passing windows are ranked by length then lower endpoint; remaining windows by descending fit p-value, length, then lower endpoint. The public tiered/first rules then choose a coherent answer.

Joint fits are unweighted. With Voigt friction, coherence requires convergence, 0<B<1, a>=0, s2>=0 and sigma_eq=a+sqrt(a*a+s2)>0. Both Voigt methods apply only the poor_fit screen; Gaussian distribution-shape checks remain not_evaluated even with check_all_windows=True. A coherent estimate failing the fit test is returned with poor_fit. Standard errors use inverse Gauss–Newton information, conditioning on the active friction boundary when a parameter is zero.

## Nonbuncher mean at the kink

`nonbuncher_mean` estimates the mean of a characteristic among nonbunchers at the kink using the observed local mixture `E[X|z] = pi(z)*mu_b + (1-pi(z))*m_NB(z)`. The buncher mean `mu_b` is fixed at the cosine-ratio `buncher_mean` estimate or a supplied value. The buncher share depends on the mass fit's Gaussian or Voigt friction law: **report the Gaussian and Voigt pair**.

Python: `specbunch.nonbuncher_mean(z, x, mass_fit, buncher_mean_fit=None, *, mu_b=None, se_mu_b=0, kink=None, weights=None, bandwidth=.05, degree=1, bin_width=.0025, threads=None)`.
R: `nonbuncher_mean(z, x, mass_fit, buncher_mean_fit=NULL, mu_b=NULL, se_mu_b=0, kink=NULL, weights=NULL, bandwidth=.05, degree=1, bin_width=.0025, threads=NULL)`.
Julia: `nonbuncher_mean(z, x, mass_fit, buncher_mean_fit=nothing; mu_b=nothing, se_mu_b=0, kink=nothing, weights=nothing, bandwidth=.05, degree=1, bin_width=.0025, threads=nothing)`.
Stata: `specbunch nonbuncher_mean z x [if] [in] [pw=weight], massfit(stored_mass) buncherfit(stored_mean) bandwidth(.05) degree(1)`; alternatively supply `mub(#) semub(#)` instead of `buncherfit()`, with advanced `binwidth(.0025) threads(#)`. Store first-stage fits with `estimates store`; omitted `massfit()` uses the current coherent mass result. The command replaces the active estimation results, and returns `e(m_NB)`, `e(se)`, `e(pi_kink)` and `e(local_fit)` with column names.

### Main options

| Option | Meaning |
|---|---|
| `z`, `x` | Equal-length running-variable and characteristic vectors; `z` is already in log units and is not logged by the package. |
| `mass_fit` | A coherent mass result fixes the atom mass and friction law; use the same population, weights and coordinate system for all stages. |
| `buncher_mean_fit` | A coherent cosine-ratio result fixes `mu_b` and supplies its SE; its kink must match the mass fit. |
| `mu_b`, `se_mu_b` | An explicit fixed buncher mean and its nonnegative SE, used instead of `buncher_mean_fit`; a zero SE treats the mean as known. |
| `kink` | The log-coordinate kink, inferred from the mass result when omitted and required to match that result when supplied. |
| `weights` | Positive observation weights, filtered jointly with finite `z` and `x`; all usable rows contribute to density normalization. |
| `bandwidth` | The triangular-kernel half-width in log units, default `.05`, chosen as a declared local sensitivity setting rather than an automatic optimum. |
| `degree` | `1` fits a local line and `2` fits a local quadratic as a sensitivity check, not a bias-corrected confidence interval. |

### Advanced (numerical settings)

`bin_width=.0025` sets bins anchored on the kink; the core smooths observed bin density minus the fixed atom density by unweighted least squares on `1, z-kink, (z-kink)^2, friction_cdf(z-kink)-.5` within `max(bandwidth,.10)`. Positive-weight bins enter the density smoother; an estimated nonbuncher density at or below `1e-12` uses that floor and is flagged. Bin-average atom probabilities divided by smoothed bin density give `pi`. The local characteristic regression uses bin weight times `max(0,1-|z-kink|/bandwidth)` and fixes the buncher mean. `threads` follows the package convention; local regressions and binning use fixed-order serial sums.

Gaussian and pure-Cauchy bin probabilities are analytic. Mixed Voigt probabilities integrate the Gaussian expectation of the Cauchy CDF over ±10 Gaussian SDs, doubling Simpson's grid from 2,000 through 64,000 intervals until successive-grid error divided by 15 is at most `1e-11`; omitted Gaussian tail probability is below `1.6e-23`. `quadrature_error` is the largest successive-grid error estimate encountered, not a rigorous total-error bound; failure to reach the criterion sets `quadrature_inconclusive`.

### Outputs and flags

| Field | Meaning |
|---|---|
| `m_NB` | Estimated nonbuncher characteristic mean at the kink, the local polynomial intercept. |
| `se` | **Standard error, not a confidence interval**: a sandwich using Kish-corrected within-bin mean variances, with buncher-mean variance propagated by the delta method assuming independence. |
| `pi_kink` | Buncher share in a bin of width `bin_width` centered at the kink, using the same smoothed density model. |
| `bandwidth`, `bin_width`, `degree`, `density_bandwidth` | The local half-width, bin width, polynomial degree and observed-density smoothing half-width actually used. |
| `mu_b`, `se_mu_b`, `B`, `sigma`, `a`, `kink` | The fixed first-stage inputs used in this calculation. |
| `n`, `n_input`, `n_eff`, `n_bins`, `threads_used` | Retained count, input count, Kish effective count, fitted local-bin count and actual runtime team size. |
| `status`, `status_code`, `refusal_reason` | `ok`/zero for an answer or `refused`/nonzero with the reason below. |
| `flags` | Semantic flag names; a finite answer is retained when a diagnostic is flagged. |
| `outside_outcome_range` | The estimated mean is outside the retained characteristic's observed range. |
| `density_floor_used` | The local fit or kink evaluation required the nonbuncher-density floor. |
| `sparse_bins` | At least one fitted bin has only one effective observation, so its within-bin variance contribution is zero. |
| `quadrature_inconclusive` | Mixed-Voigt integration exhausted the refinement budget without meeting its error criterion. |
| `local_fit` | A data frame in Python (dictionary of arrays without pandas) and R, a NamedTuple of vectors in Julia, and `e(local_fit)` in Stata. |

`local_fit` columns are `center` (bin midpoint in the original log coordinates), `observed_mean` (weighted characteristic mean), `fitted_mixture` (predicted observed mean), `nonbuncher_mean` (fitted local nonbuncher polynomial), `buncher_share` (smoothed mixture share) and `weight` (bin weight before kernel weighting).

Refusal reasons: `insufficient_data` means too few usable observations or local bins; `insufficient_local_variation` means the density or local polynomial system cannot be inverted; `numerical_failure` means the local estimate or SE is nonfinite. Invalid options are input errors. The SE conditions on the estimated mass, friction and smoothed density, omits their uncertainty and local approximation bias, and does not account for dependence between first-stage and local estimates. The target is the observed-coordinate nonbuncher CEF under the local mixture model; no deconvolution of the characteristic CEF is claimed.

## Elasticity of taxable income

The separate ETI interface reports the **Saez trapezoid** as the headline and **continuous completion** (the closed-form completion with continuous joins) as its companion. All shoulder and width calculations run in C++.

Python: `specbunch.eti.estimate(fit, z, tax_below=None, tax_above=None, *, net_below=None, net_above=None, kink=None, weights=None, integration_intervals=200, threads=None)`.
R: `eti(fit, z, tax_below=NULL, tax_above=NULL, net_below=NULL, net_above=NULL, kink=NULL, weights=NULL, integration_intervals=200, threads=NULL)`.
Julia: `SpecBunch.ETI.estimate(fit, z; tax_below=nothing, tax_above=nothing, net_below=nothing, net_above=nothing, kink=nothing, weights=nothing, integration_intervals=200, threads=nothing)`.
Stata: `specbunch eti z [if] [in] [pw=weight], massfit(stored_mass) taxbelow(#) taxabove(#)` or `netbelow(#) netabove(#)`; optional `integrationintervals(200) threads(#)`.

### Main options

`fit` is a coherent mass result supplying `B`, Gaussian SD `sigma`, Cauchy scale `a` and the selected frequency endpoints. `z` is **log income**, with the same kink coordinate and weights as that mass fit; levels must be transformed by the caller. `weights` supplies positive observation weights. `kink` defaults to the mass fit's kink and must agree if supplied. Supply exactly one complete pair: `tax_below`/`tax_above` are marginal tax rates (negative rates allow subsidies), or `net_below`/`net_above` are net-of-tax rates. Net rates must satisfy `net_below > net_above > 0`; `q=log(net_below/net_above)`.

### Advanced (numerical settings)

`integration_intervals=200` is the Simpson subinterval count below the selected window's lower endpoint (rounded up to even; range 4–4,000). The uniform CF grid continues through the selected upper endpoint, with at most 20,001 points. `threads` controls independent fixed chunks, reduced in a fixed order.

For friction transform `phi(t)=exp(-a*|t|-sigma^2*t^2/2)`, the density midpoint is `m = (integral_0^window_lo C(t)/phi(t) dt - B*window_lo)/pi`. The jump `J` is the weighted band mean of `t*S(t)/phi(t)`, with weights `phi(t)^2/t^2`. Shoulders are `left=m-J/2`, `right=m+J/2`. Saez uses `d=2*B/(left+right)`; continuous completion uses `d=3*B/(left+right+sqrt(left*right))`. Each elasticity is `e=d/q`.

### Outputs and interpretation

`saez.e` and `saez.d` are the headline elasticity and implied bunching interval width in log units; `continuous.e` and `continuous.d` are their continuous-completion counterparts. Python/Julia expose nested `.saez` and `.continuous` objects; R exposes lists; Stata returns `e(saez_e)`, `e(saez_d)`, `e(continuous_e)` and `e(continuous_d)`. Each rule also has `status`, either `ok` or a reason it is unavailable. `shoulder_left`, `shoulder_right`, `shoulder_midpoint`, and `shoulder_jump` are the inferred latent nonbuncher density functionals, with density measured per log unit. `log_net_ratio`, `net_below`, `net_above`, `window_lo`, `window_hi`, `integration_intervals`, `B`, `sigma`, `a`, `kink`, `n`, `n_input`, `n_eff` and `threads_used` record inputs and evaluation settings.

`negative_shoulder` flags a negative inferred side; Saez is retained if the midpoint is positive, while continuous completion is unavailable. `zero_shoulder` flags a boundary shoulder: the reported continuous width is the shortest minimum-energy completion, whose minimizer is not unique. Refusal reasons are `insufficient_data` (fewer than two usable rows or unusable total weight), `friction_underflow` (the fitted band cannot be safely deconvolved), and `nonpositive_density` (the midpoint is nonpositive or the shoulder calculation is nonfinite). Invalid rate pairs/options raise input errors. No SE or confidence interval is computed.

These are **finite-resolution CF functionals**: the midpoint omits the high-frequency remainder and the jump is a weighted band projection. Interpreting them as pointwise density shoulders requires local regularity and tail-approximation assumptions. Interpreting `e` structurally additionally requires the homogeneous isoelastic shift-and-bunch model and the supplied marginal incentives. Report the platform, mass-fit flags and the Gaussian/Voigt pair; neither completion rule resolves misspecified friction or first-stage mass bias.
