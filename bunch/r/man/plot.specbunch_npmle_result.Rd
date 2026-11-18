\name{plot.specbunch_npmle_result}
\alias{plot.specbunch_npmle_result}
\title{Plot fitted latent and observed distributions against the data}
\description{Base-graphics plot of an \code{npmle} result. The density view shows the data histogram (observed share divided by bin width), the fitted observed density after friction (bunchers included) and the fitted latent density with the bunching atom excluded, with a note giving the fixed mass and friction and a thin vertical line at the kink. The CDF view shows the empirical CDF, the fitted observed CDF and the fitted latent CDF, which jumps by B at the kink.}
\usage{\method{plot}{specbunch_npmle_result}(x, kind = "both", xlab = NULL, ...)}
\arguments{
\item{x}{An \code{npmle} result.}
\item{kind}{\code{"both"} (two panels side by side), \code{"pdf"} or \code{"cdf"}.}
\item{xlab}{Horizontal-axis label; defaults to the running-variable name, or \code{"z"}.}
\item{...}{Unused.}
}
\value{\code{x}, invisibly. Draw into a file device (for example \code{png()}) to save the figure.}
