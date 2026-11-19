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

transform_xscale <- function(plot) {
  if (is.null(plot) || isFALSE(plot)) return(NULL)
  if (isTRUE(plot)) return("scaled")
  if (is.character(plot) && length(plot) == 1L && !is.na(plot)) {
    k <- tolower(plot)
    if (k %in% c("none", "false", "no", "off")) return(NULL)
    if (k %in% c("scaled", "raw")) return(k)
  }
  stop('plot must be TRUE, FALSE, "scaled", "raw" or "none"', call. = FALSE)
}

#' Plot the empirical transform and the fitted model of a mass estimate
#'
#' Empirical cosine transform with a pointwise +/-2 SE band, the fitted curve
#' B*phi(t) over the whole grid and the shaded selected window. Joint fits add
#' a second panel with t*S(t), its band and the fitted jump*phi(t). Uses base
#' graphics only.
#' @param x A \code{mass} result.
#' @param frequency \code{"scaled"} (t times the pilot scale, default) or \code{"raw"} (t).
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
plot.specbunch_mass <- function(x, frequency = "scaled", ...) {
  xscale <- transform_xscale(frequency)
  if (is.null(xscale)) return(invisible(x))
  tr <- x$transform
  if (is.null(tr) || !nrow(tr)) return(invisible(x))
  col <- c(empirical = "#4d4d4d", band = "#d4d4d4", fit = "#1f3b73", window = "#f6d58e80",
           edge = "#b8860b", zero = "#999999")
  t <- if (xscale == "scaled") tr$t_scaled else tr$t
  sfx <- if (xscale == "scaled") "_scaled" else ""
  win <- c(x$raw[[paste0("window_lo", sfx)]], x$raw[[paste0("window_hi", sfx)]])
  haswin <- length(win) == 2L && all(is.finite(win))
  joint <- !is.null(tr$sine_scaled)
  fmt <- function(v) format(signif(v, 4))
  xlab <- if (xscale == "scaled") quote("Scaled frequency" ~ t ~ hat(sigma)[pilot]) else quote("Frequency" ~ t)
  method <- if (is.null(x[["method"]])) "cosine" else x[["method"]]
  friction <- if (is.null(x[["friction"]])) "gaussian" else x[["friction"]]
  head <- paste0(toupper(substring(method, 1, 1)), substring(method, 2), " fit, ", friction, " friction")
  note1 <- if (x$refused) bquote(.(paste0(head, "; refused: ", x$reason)))
    else if (friction == "voigt") bquote(.(paste0(head, ":")) ~~ hat(B) == .(fmt(x[["B"]])) * "," ~~ hat(a) == .(fmt(x[["a"]])) * "," ~~ hat(s)^2 == .(fmt(x[["s2"]])))
    else bquote(.(paste0(head, ":")) ~~ hat(B) == .(fmt(x[["B"]])) * "," ~~ hat(sigma) == .(fmt(x[["sigma"]])))
  note2 <- if (!x$refused && !is.null(x[["jump"]])) bquote(widehat(jump) == .(fmt(x[["jump"]]))) else ""
  panels <- list(list(y = tr$cosine, se = tr$cosine_se, fit = tr$cosine_fit, note = note1,
                      ylab = quote("Cosine transform" ~ hat(C)(t)), emp = quote("Empirical" ~ hat(C)(t)),
                      fl = quote("Fit" ~ hat(B) ~ hat(phi)(t)), pos = "topright"))
  if (joint) panels[[2L]] <- list(y = tr$sine_scaled, se = tr$sine_scaled_se, fit = tr$sine_fit, note = note2,
                                  ylab = quote("Scaled sine transform" ~ t ~ hat(S)(t)),
                                  emp = quote("Empirical" ~ t ~ hat(S)(t)), fl = quote("Fit" ~ widehat(jump) ~ hat(phi)(t)),
                                  pos = "topright")
  old <- graphics::par(mfrow = c(1L, length(panels)), mar = c(3.8, 3.9, 1.6, .6), mgp = c(2.4, .55, 0),
                       tcl = -.3, las = 1, cex = .85, bty = "l", bg = "white")
  on.exit(graphics::par(old))
  for (p in panels) {
    lo <- p$y - 2 * p$se; hi <- p$y + 2 * p$se
    yl <- range(c(lo, hi, p$fit, 0), finite = TRUE)
    yl <- yl + c(0, .18) * diff(yl)
    graphics::plot(NA, xlim = range(t), ylim = yl, xaxs = "i", xlab = xlab, ylab = "")
    graphics::title(ylab = p$ylab, line = 2.7)
    if (haswin) {
      graphics::rect(win[1L], yl[1L] - diff(yl), win[2L], yl[2L] + diff(yl), col = col[["window"]], border = NA)
      graphics::abline(v = win, lty = 2, lwd = .8, col = col[["edge"]])
    }
    graphics::abline(h = 0, lwd = .6, col = col[["zero"]])
    ok <- is.finite(lo) & is.finite(hi)
    graphics::polygon(c(t[ok], rev(t[ok])), c(lo[ok], rev(hi[ok])), col = col[["band"]], border = NA)
    graphics::lines(t, p$y, col = col[["empirical"]], lwd = .8)
    graphics::points(t, p$y, col = col[["empirical"]], pch = 16, cex = .25)
    if (any(is.finite(p$fit))) graphics::lines(t, p$fit, col = col[["fit"]], lwd = 1.5)
    graphics::box(bty = "l")
    graphics::mtext(as.expression(p$note), side = 3, line = .3, adj = 0, cex = .8)
    wl <- if (haswin) paste0("Selected window [", fmt(win[1L]), ", ", fmt(win[2L]), "]") else NULL
    graphics::legend(p$pos, as.expression(c(p$emp, quote(plain("") %+-% 2 ~ SE), p$fl, wl)),
                     col = c(col[["empirical"]], col[["band"]], col[["fit"]], if (haswin) col[["window"]]),
                     lty = c(1, NA, 1, if (haswin) NA), lwd = c(.8, NA, 1.5, if (haswin) NA),
                     pch = c(16, 15, NA, if (haswin) 15), pt.cex = c(.5, 1.8, NA, if (haswin) 1.8),
                     bty = "o", box.col = NA, bg = "white", cex = .75, seg.len = 2.2)
  }
  invisible(x)
}
