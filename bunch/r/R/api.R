#' @useDynLib specbunch, .registration = TRUE
#' @importFrom stats setNames
NULL

numeric_column <- function(x, name) {
  if (!is.numeric(x) || !is.null(dim(x))) stop(name, " must be a numeric vector", call. = FALSE)
  as.double(x)
}
weight_column <- function(weights, n) {
  if (is.null(weights)) return(NULL)
  weights <- numeric_column(weights, "weights")
  if (length(weights) != n) stop("weights must have the same length as the running variable", call. = FALSE)
  weights
}
scalar <- function(x, name, integer = FALSE, max = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || abs(x) > max ||
      (integer && x != trunc(x))) stop("Invalid ", name, call. = FALSE)
  as.double(x)
}
options_base <- function(seed, threads) {
  seed <- scalar(seed, "seed", TRUE, 2^53 - 1)
  if (seed < 0) stop("seed must be non-negative", call. = FALSE)
  list(seed = seed, threads = if (is.null(threads)) 0 else scalar(threads, "threads", TRUE, 2147483647))
}
result <- function(ans, type) {
  raw <- ans[[2L]]
  out <- raw
  out$status <- ans[[1L]]
  out$refused <- out$status != 0L
  out$raw <- raw
  out$kinds <- ans[[3L]]
  out$flags <- if (is.null(raw$flags_str) || !nzchar(raw$flags_str)) character() else strsplit(raw$flags_str, ",", fixed = TRUE)[[1L]]
  out$diagnostics <- raw
  out$simple_ols <- raw[startsWith(names(raw), "simple_ols.")]
  out$simple_lad <- raw[startsWith(names(raw), "simple_lad.")]
  out$tier <- if (is.null(raw$tier)) NA_real_ else raw$tier
  windows <- raw[startsWith(names(raw), "windows.")]
  names(windows) <- substring(names(windows), 9L)
  if(type == "mass") for(k in intersect(names(windows),c("cdf_inconsistent","moments_inconsistent","no_valid_distribution","distribution_check_inconclusive","residual_cdf_inconsistent","poor_fit"))) windows[[k]] <- c("not_evaluated","fail","pass")[windows[[k]]+2L]
  if(type == "buncher_mean" && !is.null(windows$poor_fit)) windows$flags <- ifelse(windows$poor_fit>0,"poor_fit", "")
  out$windows <- as.data.frame(windows, optional = TRUE)
  out$se_scale <- if(type == "mass") "log(B)" else "estimate"
  out$se <- if(type == "mass") raw$se_log_B else raw$se
  out$window <- if(type == "mass") c(raw$window_lo,raw$window_hi) else c(raw$band_lo,raw$band_hi)
  out$band <- out$window
  out$type <- type
  class(out) <- c(paste0("specbunch_",type), "specbunch_result")
  out
}
#' Core version and build information
#' @return A list with package, core, abi and build strings or numbers.
#' @export
version <- function() {
  x <- .Call(C_info)
  list(package = "0.1.0", core = x[[1L]], abi = x[[2L]], build = x[[3L]])
}
#' Display a spectral estimate
#' @param x An estimate returned by this package.
#' @param object An estimate to summarize.
#' @param ... Additional arguments (reserved).
#' @return The estimate, invisibly; summary returns a diagnostic list.
#' @details Standard errors are model based, not confidence intervals.
#' @export
print.specbunch_result <- function(x, ...) {
  target <- switch(x$type, mass = "B", buncher_mean = "mu", rd = "tau")
  cat(x$type, if (x$refused) "(refused)" else "", "\n")
  print(data.frame(estimate = x[[target]], se = x$se, status = x$status, tier = x$tier, flags = paste(x$flags, collapse = ",")), row.names = FALSE)
  if (x$type == "mass") cat("SE is on the log(B) scale.\n")
  cat("Model-based SE; not a confidence interval.\n")
  if (nzchar(x$reason)) cat(x$reason, "\n")
  invisible(x)
}
#' @rdname print.specbunch_result
#' @export
summary.specbunch_result <- function(object, ...) {
  target <- switch(object$type, mass = "B", buncher_mean = "mu", rd = "tau")
  list(estimate = data.frame(parameter = target, estimate = object[[target]], se = object$se),
       se_scale = object$se_scale, se_note = "Model-based SE; not a confidence interval",
       status = object$status, reason = object$reason, flags = object$flags,
       tier = object$tier, window = object$window, candidate_count = nrow(object$windows))
}
#' Estimate bunching mass with ranked window diagnostics
#' @export
mass <- function(z,kink=0,weights=NULL,method=c("gaussian","voigt","gaussian_jump"),
                 policy="tiered",threads=NULL,seed=1,check_all_windows=FALSE,
                 lp_iterations=5000,lp_retry_iterations=31250,npmle=FALSE) {
  zexpr <- substitute(z)
  method <- match.arg(method);policy <- match.arg(policy,c("tiered","first"))
  z <- numeric_column(z,"z");original_z <- z;kink <- scalar(kink,"kink");if(kink!=0) z <- z-kink
  opts <- options_base(seed,threads);opts$policy <- if(policy=="first") 2L else 0L
  if(!is.logical(check_all_windows)||length(check_all_windows)!=1L||is.na(check_all_windows)) stop("check_all_windows must be logical")
  opts$check_all_windows <- if(check_all_windows) 2L else 0L
  opts$lp_iterations <- scalar(lp_iterations,"lp_iterations",TRUE,2^53-1)
  opts$lp_retry_iterations <- scalar(lp_retry_iterations,"lp_retry_iterations",TRUE,2^53-1)
  if(min(opts$lp_iterations,opts$lp_retry_iterations)<=0) stop("LP iteration budgets must be positive")
  out <- result(.Call(C_estimate,match(method,c("gaussian","voigt","gaussian_jump"))-1L,z,NULL,weight_column(weights,length(z)),opts),"mass")
  out$kink <- kink
  if(!is.logical(npmle)||length(npmle)!=1L||is.na(npmle)) stop("npmle must be logical")
  if(npmle && !out$refused) {
    out$npmle <- get("npmle",mode="function")(out,original_z,weights=weights,threads=threads,plot=FALSE)
    out$npmle$xlabel <- npmle_xlabel(zexpr)
    if(interactive()) plot.specbunch_npmle_result(out$npmle)
  }
  out
}
#' Estimate buncher mean with coherent fallback and ranked bands
#' @export
buncher_mean <- function(z,x,kink=0,weights=NULL,threads=NULL,seed=1,draws=499) {
  z <- numeric_column(z,"z");x <- numeric_column(x,"x")
  if(length(x)!=length(z)) stop("x and z must have equal lengths")
  kink <- scalar(kink,"kink");if(kink!=0) z <- z-kink
  opts <- options_base(seed,threads);opts$R <- scalar(draws,"draws",TRUE,2147483647)
  out <- result(.Call(C_estimate,3L,z,x,weight_column(weights,length(z)),opts),"buncher_mean")
  out$kink <- kink;out
}

#' Shared C++ default thread count, including OMP_NUM_THREADS.
#' @export
default_threads <- function() .Call(C_default_threads)

#' Fixed-mass, fixed-friction spline postestimation
#' @export
npmle <- function(fit,z,weights=NULL,bins=NULL,threads=NULL,core_intervals=12,grid_points=1001,
                  exact_threshold=1000,max_iterations=200,lambda_min=.01,lambda_max=10000,
                  lambda_points=15,refinement_points=9,histogram_width=NULL,
                  quadrature_tolerance=1e-9,optimality_tolerance=1e-7,plot=TRUE) {
 if(!inherits(fit,"specbunch_mass")||fit$refused) stop("npmle requires a coherent mass result")
 kind <- npmle_plot_kind(plot);draw <- !is.null(kind) && (interactive() || !missing(plot))
 xlabel <- npmle_xlabel(substitute(z))
 z <- numeric_column(z,"z");weights <- weight_column(weights,length(z))
 if(!is.null(bins)) bins <- numeric_column(bins,"bins")
 opts <- list(B=fit$B,sigma=fit$sigma,a=if(fit$method=="voigt") fit$a else 0,kink=fit$kink,
              threads=if(is.null(threads)) 0 else scalar(threads,"threads",TRUE,2147483647))
 values <- as.list(environment())
 for(k in c("core_intervals","grid_points","exact_threshold","max_iterations","lambda_points","refinement_points")) opts[[k]] <- scalar(values[[k]],k,TRUE,2147483647)
 for(k in c("lambda_min","lambda_max","quadrature_tolerance","optimality_tolerance")) opts[[k]] <- scalar(values[[k]],k)
 opts$histogram_width <- if(is.null(histogram_width)) 0 else scalar(histogram_width,"histogram_width")
 ans <- .Call(C_npmle,z,weights,bins,opts);out <- ans[[2L]];out$raw <- ans[[2L]];out$status_code <- ans[[1L]];out$kinds <- ans[[3L]]
 for(prefix in c("histogram","curves","lambda_path")) {
  cols <- out$raw[startsWith(names(out$raw),paste0(prefix,"."))];names(cols) <- substring(names(cols),nchar(prefix)+2L)
  out[[prefix]] <- as.data.frame(cols,optional=TRUE)
 }
 for(k in names(out$curves)) out[[k]] <- out$curves[[k]]
 out$xlabel <- xlabel
 class(out) <- "specbunch_npmle_result"
 if(draw) plot.specbunch_npmle_result(out,kind)
 out
}
#' Display spline postestimation
#' @export
print.specbunch_npmle_result <- function(x,...) {
 cat("NPMLE:",x$status,"; lambda =",x$lambda,"; EDF =",x$edf,"\n")
 cat("The latent density jump is weakly identified.\n");invisible(x)
}
