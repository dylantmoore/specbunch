\name{nonbuncher_mean}
\alias{nonbuncher_mean}
\title{Local Nonbuncher Mean at the Kink}
\usage{
nonbuncher_mean(z, x, mass_fit, buncher_mean_fit = NULL,
  mu_b = NULL, se_mu_b = 0, kink = NULL, weights = NULL,
  bandwidth = .05, degree = 1, bin_width = .0025, threads = NULL)
}
\arguments{
\item{z, x}{Running variable in log units and characteristic, with equal lengths.}
\item{mass_fit}{Coherent mass result fixing the mass and friction.}
\item{buncher_mean_fit}{Coherent cosine-ratio buncher mean result, or NULL.}
\item{mu_b, se_mu_b}{Explicit fixed buncher mean and its SE, instead of a result.}
\item{kink}{Log-coordinate kink, inferred from and required to match the mass fit.}
\item{weights}{Positive observation weights; invalid rows are removed jointly.}
\item{bandwidth}{Triangular-kernel half-width in log units.}
\item{degree}{One for local linear, two for a local quadratic sensitivity check.}
\item{bin_width}{Advanced: width of bins anchored on the kink.}
\item{threads}{Advanced: core thread setting, NULL for the package default.}
}
\description{Fits the local observed characteristic mixture with the buncher mean held fixed. Binning, observed-density smoothing and regression run in the C++ core. Report the Gaussian and Voigt pair because mixture shares depend on the friction law.}
\details{The default local line has triangular weights and bandwidth .05. Density smoothing uses a quadratic plus the friction CDF on a window of half-width at least .10. The SE is not a confidence interval: it conditions on mass, friction and density smoothing, and propagates buncher-mean variance assuming independence. It omits local approximation bias and first-stage dependence. See the installed REFERENCE.md for all flags and numerical settings.}
\section{Advanced (numerical settings)}{\code{bin_width} and \code{threads} control bin resolution and execution.}
\value{A postestimation result with m_NB, se, pi_kink, bandwidth, flags, status and the local_fit data frame (center, observed_mean, fitted_mixture, nonbuncher_mean, buncher_share, weight). No plot is drawn.}
