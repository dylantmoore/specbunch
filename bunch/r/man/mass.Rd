\name{mass}
\alias{mass}
\title{Spectral mass estimator}
\description{Returns a coherent estimate with diagnostic flags and a ranked window table.}
\usage{mass(z, kink = 0, weights = NULL,
     method = c("gaussian", "voigt", "gaussian_jump"), policy = "tiered",
     threads = NULL, seed = 1, check_all_windows = FALSE,
     lp_iterations = 5000, lp_retry_iterations = 31250, npmle = FALSE)}
\arguments{
\item{z}{Observed running variable.}
\item{kink}{Atom location.}
\item{weights}{Optional positive observation weights.}
\item{method}{Gaussian levels, Voigt sensitivity or joint Gaussian jump fit.}
\item{policy}{Tiered selection or first coherent window.}
\item{threads}{Thread count or NULL for the shared default.}
\item{seed}{Nonnegative bootstrap seed; default 1.}
\item{check_all_windows}{Evaluate distribution checks on every coherent window.}
\item{lp_iterations}{Initial simplex iteration budget.}
\item{npmle}{Attach default spline postestimation for a coherent mass estimate.}
\item{lp_retry_iterations}{Retry simplex iteration budget.}
}
\value{An estimate list with status, reason, flags, parameter estimates, model-based standard errors, windows data frame and raw named core outputs.}
\details{See the installed public reference for definitions of every output and check. Standard errors are not confidence intervals.}
