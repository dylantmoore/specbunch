#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
root <- path.expand(if(length(args)) args[[1L]] else Sys.getenv("SPECPKG_BUILD", file.path(tempdir(),"spectral-build")))
B <- file.path(root,"r"); .libPaths(c(file.path(B,"lib"),.libPaths()))
order <- if(length(args)>1L && args[[2L]]=="rd-first") c("specrd","specbunch") else c("specbunch","specrd")
stopifnot(!nzchar(Sys.getenv("KMP_DUPLICATE_LIB_OK")))
set.seed(12); n <- 256L
z <- rnorm(n); d <- as.double(z>0); y <- d + rnorm(n)
first <- list()
for(t in c(1L,2L,4L,8L)) for(pkg in order) {
  library(pkg,character.only=TRUE)
  ans <- if(pkg=="specrd") specrd::rd(z,d,y,threads=t) else specbunch::buncher_mean(z,y,threads=t,R=19)
  if(is.null(first[[pkg]])) first[[pkg]] <- ans$raw else stopifnot(identical(first[[pkg]],ans$raw))
  stopifnot(!nzchar(Sys.getenv("KMP_DUPLICATE_LIB_OK")))
}
cat(paste(order,collapse=" -> "),": interleaved calls 1/2/4/8 PASS\n")
