\name{buncher_mean}
\alias{buncher_mean}
\title{Spectral buncher mean estimator}
\description{Returns a coherent estimate with diagnostic flags and a ranked window table.}
\usage{buncher_mean(z, x, kink = 0, weights = NULL, threads = NULL,
             seed = 1, draws = 499)}
\arguments{
\item{z}{Observed running variable.}
\item{x}{Outcome variable.}
\item{kink}{Atom location.}
\item{weights}{Optional positive observation weights.}
\item{threads}{Thread count or NULL for the shared default.}
\item{seed}{Nonnegative bootstrap seed; default 1.}
\item{draws}{Fit-test bootstrap count, at least 19.}
}
\value{An estimate list with status, reason, flags, parameter estimates, model-based standard errors, windows data frame and raw named core outputs.}
\details{See the installed public reference for definitions of every output and check. Standard errors are not confidence intervals.}
