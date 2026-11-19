{smcl}
{title:specbunch}

{p}Bunching mass and buncher mean with ranked sensitivity tables.{p_end}
{cmd:specbunch mass z [if] [in] [pw=weight], method(cosine) friction(gaussian) policy(tiered) [plot(scaled|raw) noplot name(name[, replace]) saving(filename[, replace])]}
{cmd:specbunch buncher z x [if] [in] [pw=weight], draws(499)}
{cmd:specbunch windows}

{p}Methods: cosine (default), joint. Friction: gaussian (default), voigt. Rules: tiered (default), first. Options: kink(real), threads(integer), seed(integer; default 1), checkallwindows, lpiterations(integer; default 5000), lpretryiterations(integer; default 31250). joint requires unweighted data. Mean draws must be at least 19.{p_end}
{p}e(windows) contains ranked sensitivity columns. Distribution checks use 1=pass, 0=fail, -1=not evaluated. Coherent estimates are returned with flags; refusal means insufficient_data, no_frequency_band, or no_coherent_window.{p_end}
{p}Flags: cdf_inconsistent (DKW), moments_inconsistent (harmonic bounds), no_valid_distribution (LP exclusion), distribution_check_inconclusive (LP budget), residual_cdf_inconsistent (probit screen), poor_fit (p below .10).{p_end}
{p}See docs/REFERENCE.md for the full option and result reference. Standard errors are model based, not confidence intervals.{p_end}

{title:Spline postestimation}
{cmd:specbunch npmle [, bins(matrix_name) threads(#) coreintervals(12) gridpoints(1001) exactthreshold(1000) maxiterations(200) lambdamin(.01) lambdamax(10000) lambdapoints(15) refinepoints(9) histogramwidth(0) quadraturetolerance(1e-9) optimalitytolerance(1e-7) plot(both|pdf|cdf) noplot name(name[, replace]) saving(filename[, replace])]}
{p}B and friction are fixed at the preceding mass estimate. The full retained sample contributes. Default likelihood bins have width max(sigma,a)/20; output bins do not restrict the fit. BIC selects the spline penalty with interior refinement; e(npmle_lambda_at_grid_edge) identifies a capped endpoint minimum.{p_end}
{p}e(npmle_curves) columns are grid, latent_pdf, latent_cdf, observed_pdf, observed_cdf and empirical_cdf. The latent PDF excludes the atom; the latent CDF includes its jump. e(npmle_histogram) reports bin_lo, bin_hi, observed_share, fitted_share, nonbuncher_share, buncher_share and latent_counterfactual_density.{p_end}
{p}By default {cmd:specbunch npmle} draws a graph. The density view shows the data histogram (observed share divided by bin width), the fitted observed density after friction (bunchers included) and the fitted latent density excluding the bunching atom, with a note giving B and friction and a dotted line at the kink. The CDF view shows the empirical, fitted observed and fitted latent CDF; the latent CDF jumps by B at the kink. {opt plot(both)} (the default) combines the two views with {cmd:graph combine}; {opt plot(pdf)} or {opt plot(cdf)} draws one; {opt noplot} skips drawing. {opt name()} and {opt saving()} are passed to the final graph. The horizontal range shows the 0.1% to 99.9% empirical quantiles; this affects the display only. The data and e(sample) are unchanged.{p_end}
{p}Scalars use the e(npmle_) prefix: lambda, edf, loglik, bic, latent_left, latent_right, latent_jump and numerical diagnostics. e(npmle_status) reports convergence. The latent density jump is weakly identified. See the full reference for every field, option and error convention.{p_end}

{title:NPMLE main options}
{p}{opt bins()} sets output edges; weights are inherited from the mass sample. {opt plot(both|pdf|cdf)}, {opt noplot}, {opt name()} and {opt saving()} control drawing.{p_end}

{title:Advanced (solver settings)}
{p}{opt threads()}, {opt coreintervals()}, {opt gridpoints()}, {opt exactthreshold()}, {opt maxiterations()}, {opt lambdamin()}, {opt lambdamax()}, {opt lambdapoints()}, {opt refinepoints()}, {opt histogramwidth()}, {opt quadraturetolerance()} and {opt optimalitytolerance()} control numerical evaluation. The initial lambda grid extends automatically until BIC has an interior minimum or a search cap is reached; see REFERENCE.md for defaults and meanings.{p_end}

{p}Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.{p_end}

{p}e(transform) contains t, t_scaled, cosine, cosine_se and cosine_fit; joint fits add sine_scaled, sine_scaled_se and sine_fit. e(window_lo), e(window_hi) and e(n_points) describe the selected fitting window. See REFERENCE.md for the pointwise variance formula.{p_end}
{p}By default {cmd:specbunch mass} draws the transform fit: the empirical cosine transform with a pointwise ±2 SE band, the fitted curve over the whole grid and the selected window shaded between dashed lines, titled with method, friction and the estimates. Joint fits add a panel with t times the sine transform and its fit, combined with {cmd:graph combine}. {opt plot(scaled)} (the default) uses frequency times the pilot scale on the horizontal axis; {opt plot(raw)} uses frequency t; {opt noplot} skips drawing. {opt name()} and {opt saving()} are passed to the final graph; give the mass and NPMLE graphs different names to keep both. Drawing uses preserve/restore.{p_end}

{p}All four method/friction combinations are supported. Joint fits are unweighted. Voigt fits use coherence and poor_fit only; Gaussian distribution-shape checks are not applied.{p_end}

{title:Nonbuncher mean at the kink}
{p 4 4 2}{cmd:specbunch nonbuncher_mean z x} {it:[if] [in] [pw=weight]},
{cmd:massfit(stored_mass)} {cmd:buncherfit(stored_mean)}
[{cmd:bandwidth(.05) degree(1) binwidth(.0025) threads(0)}].
Store the mass and buncher-mean estimates with {cmd:estimates store}.
Alternatively use {cmd:mub(#) semub(#)} instead of {cmd:buncherfit()}.
Omitted {cmd:massfit()} uses the current coherent mass result.
Both first stages must use the same kink, weights and population.
{p 4 4 2}The C++ core bins raw log-coordinate z and fits the observed mixture
with the buncher mean fixed. Triangular local-linear fitting at bandwidth .05
is the default; {cmd:degree(2)} is a quadratic sensitivity check.
Report the Gaussian and Voigt pair: mixture shares depend on friction.
{p 4 4 2}Returns {cmd:e(m_NB)}, {cmd:e(se)}, {cmd:e(pi_kink)}, settings and flags,
and {cmd:e(local_fit)} with center, observed_mean, fitted_mixture,
nonbuncher_mean, buncher_share and weight columns. The command replaces e().
SE is not a confidence interval: it conditions on mass, friction and smoothed
density and propagates buncher-mean variance assuming independence.
{p 4 4 2}Flags: outside_outcome_range (estimate outside observed X range),
density_floor_used (nonbuncher density required a positive floor), sparse_bins
(a fitted bin has only one effective record), quadrature_inconclusive
(Voigt integration did not meet its refinement criterion).
Refusals: insufficient_data, insufficient_local_variation, numerical_failure.
{p 4 4 2}{bf:Advanced numerical settings:} {cmd:binwidth()} sets the bins;
{cmd:threads()} follows the core convention. See REFERENCE.md for density
smoothing and conditional-SE assumptions. No plot is drawn.

{title:Elasticity of taxable income}
{p 4 4 2}{cmd:specbunch eti z} {it:[if] [in] [pw=weight]},
{cmd:massfit(stored_mass) taxbelow(#) taxabove(#)}; alternatively supply
{cmd:netbelow(#) netabove(#)}. Supply exactly one rate pair with
netbelow greater than netabove greater than zero. Negative tax rates permit
subsidies. Input z must already be in log units, with the mass fit's kink.
{p 4 4 2}The Saez trapezoid is the headline; continuous completion is the
companion. The core computes finite-resolution CF shoulders and widths,
then e=d/log(netbelow/netabove). Returns {cmd:e(saez_e)}, {cmd:e(saez_d)},
{cmd:e(continuous_e)}, {cmd:e(continuous_d)}, {cmd:e(shoulder_left)},
{cmd:e(shoulder_right)}, midpoint, jump and flags. No SE or CI is computed.
{p 4 4 2}negative_shoulder flags a negative side and makes continuous
completion unavailable, retaining Saez when its midpoint is positive.
zero_shoulder flags a boundary where the shortest continuous completion is
reported. Refusals: insufficient_data, friction_underflow, nonpositive_density.
{p 4 4 2}{bf:Advanced numerical settings:} {cmd:integrationintervals(200)}
sets Simpson intervals below the mass window; {cmd:threads(0)} controls
fixed-chunk execution. Structural interpretation requires homogeneous
isoelastic shifting and the supplied marginal incentives; the shoulders are
finite-resolution functionals. Report Gaussian/Voigt results and mass-fit flags.

{p}NPMLE accepts {opt basis(cubic|linear)} and optional taxbelow()/taxabove() or netbelow()/netabove(). When rates are supplied, ETI always reports the transform and NPMLE shoulder pair. ETI also automatically uses an attached NPMLE, or fits one with {opt npmle}. Additional scalars use e(npmle_saez_e), e(npmle_continuous_e), and e(npmle_shoulder_left/right). Grid-edge, shoulders_weak, and shoulders_inconclusive diagnostics are reported; see REFERENCE.md for the BIC-nearby sensitivity definition. No confidence intervals are computed.{p_end}

{title:Licence and citing}
{pstd}Own code: MIT, Copyright (c) 2026 Dylan Moore. See specbunch_THIRD_PARTY_NOTICES.txt and supplied licences.
{pstd}Please cite Moore, Dylan (2026). An Application of Euler's Formula to Bunching Methodology: Addressing Frictions via Spectral Deconvolution. Working paper. This title is a placeholder to be updated.
