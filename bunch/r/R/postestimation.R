post_mass_options <- function(fit,kink,threads) {
 if(!inherits(fit,"specbunch_mass") || fit$refused) stop("a coherent mass result is required")
 if(is.null(kink)) kink <- fit$kink
 kink <- scalar(kink,"kink")
 if(kink != fit$kink) stop("kink must match the mass fit")
 list(B=fit$B,sigma=fit$sigma,a=if(fit$friction=="voigt") fit$a else 0,kink=kink,
      threads=if(is.null(threads)) 0 else scalar(threads,"threads",TRUE,2147483647))
}
post_result <- function(ans) {
 out <- ans[[2L]];out$raw <- ans[[2L]];out$kinds <- ans[[3L]];out$status_code <- ans[[1L]];out$refused <- out$status_code!=0L
 out$flags <- if(nzchar(out$raw$flags)) strsplit(out$raw$flags," ",fixed=TRUE)[[1L]] else character()
 for(prefix in c("local_fit","saez","continuous")) {
  cols <- out$raw[startsWith(names(out$raw),paste0(prefix,"."))];names(cols) <- substring(names(cols),nchar(prefix)+2L)
  out[[prefix]] <- if(prefix=="local_fit") as.data.frame(cols,optional=TRUE) else cols
 }
 class(out) <- "specbunch_postestimation";out
}
#' Local nonbuncher mean with a fixed buncher mean
#' @export
nonbuncher_mean <- function(z,x,mass_fit,buncher_mean_fit=NULL,mu_b=NULL,se_mu_b=0,
                           kink=NULL,weights=NULL,bandwidth=.05,degree=1,bin_width=.0025,threads=NULL) {
 opts <- post_mass_options(mass_fit,kink,threads)
 if(!is.null(buncher_mean_fit)) {
  if(!is.null(mu_b)) stop("supply buncher_mean_fit or mu_b, not both")
  if(!inherits(buncher_mean_fit,"specbunch_buncher_mean") || buncher_mean_fit$refused) stop("a coherent buncher_mean result is required")
  if(buncher_mean_fit$kink!=opts$kink) stop("buncher mean and mass fit must have the same kink")
  mu_b <- buncher_mean_fit$mu;se_mu_b <- buncher_mean_fit$se
 }
 if(is.null(mu_b)) stop("buncher_mean_fit or mu_b is required")
 z <- numeric_column(z,"z");x <- numeric_column(x,"x");if(length(z)!=length(x)) stop("x and z must have equal lengths")
 opts$mu_b <- scalar(mu_b,"mu_b");opts$se_mu_b <- scalar(se_mu_b,"se_mu_b")
 opts$bandwidth <- scalar(bandwidth,"bandwidth");opts$degree <- scalar(degree,"degree",TRUE,2147483647);opts$bin_width <- scalar(bin_width,"bin_width")
 post_result(.Call(C_post,0L,z,x,weight_column(weights,length(z)),opts))
}
#' Elasticity of taxable income by shift-and-complete
#' @export
eti <- function(fit,z,tax_below=NULL,tax_above=NULL,net_below=NULL,net_above=NULL,
                kink=NULL,weights=NULL,integration_intervals=200,threads=NULL) {
 opts <- post_mass_options(fit,kink,threads)
 if(is.null(net_below)&&is.null(net_above)) {
  if(is.null(tax_below)||is.null(tax_above)) stop("supply both tax rates or both net-of-tax rates")
  net_below <- 1-scalar(tax_below,"tax_below");net_above <- 1-scalar(tax_above,"tax_above")
 } else if(!is.null(tax_below)||!is.null(tax_above)||is.null(net_below)||is.null(net_above)) stop("supply exactly one complete rate pair")
 opts$net_below <- scalar(net_below,"net_below");opts$net_above <- scalar(net_above,"net_above")
 opts$window_lo <- fit$window_lo;opts$window_hi <- fit$window_hi
 opts$integration_intervals <- scalar(integration_intervals,"integration_intervals",TRUE,2147483647)
 z <- numeric_column(z,"z");post_result(.Call(C_post,1L,z,NULL,weight_column(weights,length(z)),opts))
}
#' Display postestimation results
#' @export
print.specbunch_postestimation <- function(x,...) {
 cat(x$estimator,": ",x$status,"\n",sep="")
 if(!is.null(x$m_NB)) {cat("m_NB =",x$m_NB,"; SE =",x$se,"\n");cat(x$se_note,"\n")} else {print(x$saez);print(x$continuous)}
 if(length(x$flags)) cat("Flags:",paste(x$flags,collapse=", "),"\n")
 invisible(x)
}
