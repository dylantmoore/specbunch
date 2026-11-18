\name{mass}
\alias{mass}
\title{Spectral mass estimator}
\description{Returns a coherent estimate with diagnostic flags and a ranked window table.}
\usage{mass(z, kink = 0, weights = NULL,
     method = c("cosine", "joint"), friction = c("gaussian", "voigt"), policy = "tiered",
     threads = NULL, seed = 1, check_all_windows = FALSE,
     lp_iterations = 5000, lp_retry_iterations = 31250, npmle = FALSE)}
\arguments{
\item{z}{Observed running variable.}
\item{kink}{Atom location.}
\item{weights}{Optional positive observation weights.}
\item{method}{Cosine-transform fit, or joint cosine and sine fit allowing a density jump.}
\item{friction}{Gaussian or symmetric Voigt friction, for either method.}
\item{policy}{Tiered selection or first coherent window.}
\item{threads}{Thread count or NULL for the shared default.}
\item{seed}{Nonnegative bootstrap seed; default 1.}
\item{check_all_windows}{Evaluate distribution checks on every coherent window.}
\item{lp_iterations}{Initial simplex iteration budget.}
\item{npmle}{Attach default spline postestimation for a coherent mass estimate.}
\item{lp_retry_iterations}{Retry simplex iteration budget.}
}
\value{An estimate list with status, reason, flags, parameter estimates, model-based standard errors, windows data frame and raw named core outputs.}
\details{All four method/friction combinations are supported. Joint fits are unweighted. Voigt fits use coherence and poor_fit only; Gaussian distribution-shape checks are not applied.
See the installed public reference for definitions of every output and check. Standard errors are not confidence intervals.

Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.

The transform table contains t, t_scaled, cosine, cosine_se and cosine_fit; joint fits also return sine_scaled, sine_scaled_se and sine_fit. The selected window endpoints and n_points are separate scalars. See the installed REFERENCE.md for column definitions.}
