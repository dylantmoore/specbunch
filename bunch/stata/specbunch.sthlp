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
