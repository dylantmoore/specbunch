npmle_plot_kind <- function(plot) {
  if (is.null(plot) || isFALSE(plot)) return(NULL)
  if (isTRUE(plot)) return("both")
  if (is.character(plot) && length(plot) == 1L && !is.na(plot)) {
    k <- tolower(plot)
    if (k %in% c("none", "false", "no", "off")) return(NULL)
    if (k %in% c("both", "pdf", "cdf")) return(k)
  }
  stop('plot must be TRUE, FALSE, "both", "pdf", "cdf" or "none"', call. = FALSE)
}

npmle_xlabel <- function(expr) if (is.symbol(expr)) as.character(expr) else "z"

npmle_view_range <- function(grid, ecdf, kink) {
  at <- function(p) grid[min(which(ecdf >= p)[1L], length(grid), na.rm = TRUE)]
  lo <- min(at(.001), kink); hi <- max(at(.999), kink)
  pad <- .04 * if (hi > lo) hi - lo else 1
  c(max(lo - pad, grid[1L]), min(hi + pad, grid[length(grid)]))
}

#' Plot fitted latent and observed distributions against the data
#'
#' Density view: data histogram (observed share / bin width), fitted observed
#' density (after friction, bunchers included) and fitted latent density
#' excluding the bunching atom. CDF view: empirical, fitted observed and fitted
#' latent CDF, which jumps by B at the kink. Uses base graphics only.
#' @param x An \code{npmle} result.
#' @param kind \code{"both"} (default), \code{"pdf"} or \code{"cdf"}.
#' @param xlab Axis label; defaults to the running-variable name or \code{"z"}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
plot.specbunch_npmle_result <- function(x, kind = "both", xlab = NULL, ...) {
  kind <- npmle_plot_kind(kind)
  if (is.null(kind)) return(invisible(x))
  col <- c(data = "#b3b3b3", observed = "#1f3b73", latent = "#d9822b", kink = "#555555")
  cv <- x$curves; h <- x$histogram; kink <- x[["kink"]]
  if (is.null(xlab)) xlab <- if (is.null(x[["xlabel"]])) "z" else x[["xlabel"]]
  xlim <- npmle_view_range(cv$grid, cv$empirical_cdf, kink)
  panels <- if (kind == "both") c("pdf", "cdf") else kind
  old <- graphics::par(mfrow = c(1L, length(panels)), mar = c(3.6, 3.8, 1, .6), mgp = c(2.3, .55, 0),
                       tcl = -.3, las = 1, cex = .85, bty = "l", bg = "white")
  on.exit(graphics::par(old))
  fmt <- function(v) formatC(v, digits = 4, format = "g")
  for (panel in panels) {
    if (panel == "pdf") {
      dens <- h$observed_share / (h$bin_hi - h$bin_lo)
      inview <- h$bin_hi >= xlim[1L] & h$bin_lo <= xlim[2L]
      top <- max(c(dens[inview], cv$observed_pdf, cv$latent_pdf), na.rm = TRUE) * 1.12
      graphics::plot(NA, xlim = xlim, ylim = c(0, top), xaxs = "i", yaxs = "i", xlab = xlab, ylab = "Density")
      graphics::rect(h$bin_lo, 0, h$bin_hi, dens, col = col[["data"]], border = NA)
      graphics::abline(v = kink, lty = 3, lwd = .8, col = col[["kink"]])
      graphics::lines(cv$grid, cv$observed_pdf, col = col[["observed"]], lwd = 1.3)
      graphics::lines(cv$grid, cv$latent_pdf, col = col[["latent"]], lwd = 1.3, lty = 2)
      note <- c(bquote(hat(B) == .(fmt(x[["B"]]))), bquote(hat(sigma) == .(fmt(x[["sigma"]]))))
      if (isTRUE(x[["a"]] > 0)) note <- c(note, bquote(hat(a) == .(fmt(x[["a"]]))))
      if (!identical(x$status, "converged")) note <- c(note, paste("status:", x$status))
      graphics::legend("topleft", legend = as.expression(note), bty = "n", cex = .9, inset = c(-.02, 0))
      graphics::legend("topright", c("Data (histogram)", "Observed fit", "Latent fit, excl. bunchers", "Kink"),
                       col = col[c("data", "observed", "latent", "kink")], lty = c(NA, 1, 2, 3),
                       lwd = c(NA, 1.3, 1.3, .8), pch = c(15, NA, NA, NA), pt.cex = 1.6, bty = "o", box.col = NA, bg = "white", cex = .75,
                       seg.len = 2.2)
    } else {
      graphics::plot(NA, xlim = xlim, ylim = c(0, 1.02), xaxs = "i", yaxs = "i", xlab = xlab,
                     ylab = "Cumulative probability")
      graphics::abline(v = kink, lty = 3, lwd = .8, col = col[["kink"]])
      graphics::lines(cv$grid, cv$empirical_cdf, type = "s", col = col[["data"]], lwd = 2.6)
      graphics::lines(cv$grid, cv$observed_cdf, col = col[["observed"]], lwd = 1.3)
      graphics::lines(cv$grid, cv$latent_cdf, col = col[["latent"]], lwd = 1.3, lty = 2)
      graphics::legend("topleft", as.expression(c("Empirical CDF", "Observed fit",
                                                  bquote("Latent fit (jump" ~ hat(B) ~ "at kink)"), "Kink")),
                       col = col[c("data", "observed", "latent", "kink")], lty = c(1, 1, 2, 3),
                       lwd = c(2.6, 1.3, 1.3, .8), bty = "n", cex = .75, seg.len = 2.2)
    }
  }
  invisible(x)
}
