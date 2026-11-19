\name{npmle}
\alias{npmle}
\title{Fixed-parameter spline postestimation}
\description{Convex penalised likelihood with nonnegative cubic splines on each side of the kink, fixed mass and friction, and BIC-selected smoothing.}
\usage{npmle(fit, z, weights = NULL, bins = NULL, threads = NULL,
      core_intervals = 12, grid_points = 1001, exact_threshold = 1000,
      max_iterations = 200, lambda_min = .01, lambda_max = 10000,
      lambda_points = 15, refinement_points = 9, histogram_width = NULL,
      quadrature_tolerance = 1e-9, optimality_tolerance = 1e-7, plot = TRUE, basis = "cubic", tax_below = NULL, tax_above = NULL,
      net_below = NULL, net_above = NULL)}
\arguments{
\item{basis}{Clamped spline basis: "cubic" (default) or "linear".}
\item{tax_below, tax_above, net_below, net_above}{Optional complete tax or net-rate pair; automatically attaches both ETI shoulder calculations.}
\item{fit}{Coherent mass result supplying fixed mass, friction and kink.}
\item{z}{Full running-variable sample.}
\item{weights}{Optional positive observation weights.}
\item{bins}{Optional finite increasing output-bin edges; does not restrict the likelihood sample.}
\item{threads}{Core thread count.}
\item{core_intervals}{Central spline knot intervals per side.}
\item{grid_points}{Central curve-grid point count, with additional full-range tail points.}
\item{exact_threshold}{Largest sample evaluated observation by observation.}
\item{max_iterations}{Newton iteration budget.}
\item{lambda_min, lambda_max}{Initial penalty-grid endpoints; extended automatically until an interior BIC minimum or the search caps.}
\item{lambda_points}{Coarse penalty-grid count.}
\item{refinement_points}{Interior refinement count.}
\item{histogram_width}{Likelihood-bin width; NULL uses max(sigma,a)/20.}
\item{quadrature_tolerance}{Target scaled quadrature error.}
\item{optimality_tolerance}{Normalized simplex KKT convergence tolerance.}
\item{plot}{\code{TRUE} (default, same as \code{"both"}), \code{"pdf"}, \code{"cdf"}, or \code{FALSE}/\code{"none"}. The default draws with base graphics only in interactive sessions; an explicit value also draws in batch. See \code{\link{plot.specbunch_npmle_result}}.}
}
\value{A list containing status, penalty, EDF, log-likelihood, BIC, latent density limits, histogram and curve data frames, and raw core outputs.}
\details{The latent jump is weakly identified. Curves contain grid, latent_pdf, latent_cdf, observed_pdf, observed_cdf and empirical_cdf. The latent PDF excludes the atom; its CDF includes it. By default an interactive call draws the fitted distributions against the data; \code{plot(result)} redraws them. The installed REFERENCE.md describes every output and numerical convention.}

\section{Main options}{\code{weights}, \code{bins} and \code{plot} control the sample weights, output histogram and drawing.}
\section{Advanced (solver settings)}{\code{threads}, \code{core_intervals}, \code{grid_points}, \code{exact_threshold}, \code{max_iterations}, \code{lambda_min}, \code{lambda_max}, \code{lambda_points}, \code{refinement_points}, \code{histogram_width}, \code{quadrature_tolerance} and \code{optimality_tolerance} control numerical evaluation; see their argument descriptions and the installed REFERENCE.md.}
