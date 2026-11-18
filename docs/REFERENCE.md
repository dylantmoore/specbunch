# Public reference

`specbunch` estimates a bunching mass and a buncher outcome mean using frequency-domain fits. Python, R, Julia and Stata call the same deterministic C++ core. NPMLE postestimation provides default density and CDF plots through the language wrappers.

## Mass

The headline estimator is weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits the sine transform. The companions are log-linear OLS and LAD on the reference band.

Python: `mass(z, kink=0, weights=None, method="cosine", friction="gaussian", policy="tiered", threads=None, seed=1, check_all_windows=False, lp_iterations=5000, lp_retry_iterations=31250, npmle=False)`.
R uses the same names with `NULL` and `FALSE`; Julia uses keyword arguments and method/friction/policy symbols. Stata uses `specbunch mass z [if] [in] [pw=weight], method(cosine) friction(gaussian) policy(tiered)`.

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

Every mass result includes a numeric `transform` table (Python pandas DataFrame when available, otherwise a dict of NumPy arrays; R data.frame; Julia NamedTuple of vectors; Stata e(transform)). No transform plotting code is included.
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

## Joint transform fit with Voigt friction

All four method/friction combinations are supported. For symmetric Voigt friction, `phi(t) = exp(-a*abs(t) - s2*t*t/2)`; the joint model fits `B*phi(t)` to the empirical cosine transform and `jump*phi(t)/t` to the empirical sine transform. Convolution by symmetric friction multiplies both Fourier components by the same real characteristic function, so the jump term carries the same phi. The covariance weights are fixed from the empirical cosine/sine covariance at each frequency.

The joint Voigt solver profiles the two linear coefficients (B and the scaled jump) exactly, then uses constrained Gauss–Newton steps with line search for nonnegative a and s2, checking both friction boundaries. Ranked windows use the joint Gaussian frequency-grid and visibility-extension design: at least 11 frequencies, a pilot-scaled spacing of 0.05 from 0.2 to an initial upper bound of 4, extended up to 8 when signal remains visible. Its fit test uses 499 deterministic sign-multiplier draws, refitting the joint model and comparing the maximum standardized cosine/sine residual. Passing windows are ranked by length then lower endpoint; remaining windows by descending fit p-value, length, then lower endpoint. The public tiered/first rules then choose a coherent answer.

Joint fits are unweighted. With Voigt friction, coherence requires convergence, 0<B<1, a>=0, s2>=0 and sigma_eq=a+sqrt(a*a+s2)>0. Both Voigt methods apply only the poor_fit screen; Gaussian distribution-shape checks remain not_evaluated even with check_all_windows=True. A coherent estimate failing the fit test is returned with poor_fit. Standard errors use inverse Gauss–Newton information, conditioning on the active friction boundary when a parameter is zero.
