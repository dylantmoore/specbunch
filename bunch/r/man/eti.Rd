\name{eti}
\alias{eti}
\title{Elasticity of Taxable Income by Shift and Completion}
\usage{
eti(fit, z, tax_below = NULL, tax_above = NULL,
  net_below = NULL, net_above = NULL, kink = NULL, weights = NULL,
  integration_intervals = 200, threads = NULL)
}
\arguments{
\item{fit}{Coherent mass result supplying mass, friction and frequency window.}
\item{z}{Running variable already expressed in log units.}
\item{tax_below, tax_above}{Marginal tax rates on the two sides; subsidies may be negative.}
\item{net_below, net_above}{Alternative net-of-tax rates, with net_below greater than net_above greater than zero; supply exactly one complete rate pair.}
\item{kink}{Log-coordinate kink, inferred from and required to match the fit.}
\item{weights}{Positive observation weights.}
\item{integration_intervals}{Advanced: Simpson interval count, rounded up to even.}
\item{threads}{Advanced: core thread setting.}
}
\description{Reports the Saez trapezoid as the headline and continuous completion as its companion, using characteristic-function density shoulders computed in C++.}
\details{The midpoint is a finite-frequency deconvolution integral; the jump is a weighted sine-band projection. These are finite-resolution functionals, not unrestricted pointwise densities. Structural interpretation requires homogeneous isoelastic responses and the specified marginal incentives. Negative shoulders flag the answer and make continuous completion unavailable; Saez remains when the midpoint is positive. No SE or confidence interval is computed. Report the Gaussian/Voigt pair and first-stage flags.}
\section{Advanced (numerical settings)}{\code{integration_intervals} defaults to 200 below the selected lower frequency; \code{threads} controls fixed-chunk execution.}
\value{A postestimation result with nested saez and continuous lists (e, d, status), shoulder_left, shoulder_right, shoulder_midpoint, shoulder_jump, flags and input diagnostics.}
