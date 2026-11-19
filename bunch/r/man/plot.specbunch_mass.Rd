\name{plot.specbunch_mass}
\alias{plot.specbunch_mass}
\title{Plot the empirical transform and the fitted model of a mass estimate}
\description{Empirical cosine transform with a pointwise +/-2 SE band, the fitted curve B*phi(t) over the whole grid and the shaded selected window, with a note giving method, friction and estimates. Joint fits add a second panel with t*S(t), its band and the fitted jump*phi(t). Uses base graphics only.}
\usage{\method{plot}{specbunch_mass}(x, frequency = "scaled", ...)}
\arguments{
\item{x}{A \code{mass} result.}
\item{frequency}{\code{"scaled"} (t times the pilot scale, default) or \code{"raw"} (t).}
\item{...}{Unused.}
}
\value{\code{x}, invisibly.}
