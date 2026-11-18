library(specbunch)
serial_build <- grepl("openmp=off", version()$build, fixed=TRUE)
used <- function(n) if(serial_build) 1 else n
set.seed(42)
n <- 256L
z <- ifelse(runif(n) < .2, 0, rnorm(n)) + .04*rnorm(n)
x <- 2 + .2*rnorm(n)
for(method in c("gaussian","voigt","gaussian_jump")) {
  a <- mass(z,method=method,threads=1)
  b <- mass(c(z,NA_real_,Inf),method=method,threads=2)
  a$raw$n_input <- b$raw$n_input
stopifnot(a$threads_used == 1, b$threads_used %in% used(c(2,8)))
a$raw$threads_used <- b$raw$threads_used
  stopifnot(identical(a$raw,b$raw),is.data.frame(a$windows),inherits(a,"specbunch_result"),identical(a$se,a$raw$se_log_B))
  print(a); stopifnot(is.list(summary(a)))
}
a <- buncher_mean(z,x,threads=1,draws=19)
b <- buncher_mean(c(z,NA_real_),c(x,0),threads=8,draws=19)
a$raw$n_input <- b$raw$n_input
stopifnot(a$threads_used == 1, b$threads_used %in% used(c(2,8)))
a$raw$threads_used <- b$raw$threads_used
stopifnot(identical(a$raw,b$raw),is.list(version()))
errors <- list(quote(mass(z,seed=-1)),quote(mass(z,weights=1)),
               quote(mass(z,method="gaussian_jump",weights=rep(1,n))),quote(buncher_mean(z,x[-1])))
for(e in errors) stopifnot(inherits(try(eval(e),silent=TRUE),"try-error"))
stopifnot(identical(mass(z+1,kink=1)$raw$n,as.double(n)))

old_omp <- Sys.getenv("OMP_NUM_THREADS", unset=NA_character_)
Sys.unsetenv("OMP_NUM_THREADS")
stopifnot(buncher_mean(z,x,draws=19)$threads_used == used(default_threads()))
for(t in c(1,3,16)) stopifnot(buncher_mean(z,x,draws=19,threads=t)$threads_used == used(t))
Sys.setenv(OMP_NUM_THREADS="3")
stopifnot(default_threads()==3L, buncher_mean(z,x,draws=19)$threads_used==used(3), buncher_mean(z,x,draws=19,threads=16)$threads_used==used(16))
if(is.na(old_omp)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS=old_omp)

stopifnot(mass(1:10)$reason == "insufficient_data")

set.seed(37)
nz <- ifelse(runif(1500)<.2,0,rnorm(1500))+.04*rnorm(1500)
nf <- mass(nz,threads=1,npmle=TRUE)
stopifnot(!nf$refused,inherits(nf$npmle,"specbunch_npmle_result"))
nr <- npmle(nf,nz,threads=1,grid_points=101,lambda_points=5,refinement_points=3)
stopifnot(nr$status=="converged",is.data.frame(nr$curves),is.data.frame(nr$histogram),all(nr$coefficients>=0),abs(tail(nr$empirical_cdf,1)-1)<1e-12)
stopifnot(all(diff(nr$observed_cdf)>=-1e-12),all(diff(nr$latent_cdf)>=-1e-12))
stopifnot(identical(nr$histogram$fitted_share,nr$histogram$nonbuncher_share+nr$histogram$buncher_share))
print(nr)

pf <- tempfile(fileext=".png");grDevices::png(pf,width=1500,height=585,res=150)
dev_before <- grDevices::dev.cur()
pr <- npmle(nf,nz,threads=1,grid_points=101,lambda_points=3,refinement_points=3,plot="both")
stopifnot(identical(pr$xlabel,"nz"),identical(plot(pr,kind="cdf"),pr))
grDevices::dev.off();stopifnot(file.exists(pf),file.info(pf)$size>5000)
n_dev <- length(grDevices::dev.list())
quiet <- npmle(nf,nz,threads=1,grid_points=101,lambda_points=3,refinement_points=3,plot=FALSE)
stopifnot(length(grDevices::dev.list())==n_dev)
stopifnot(inherits(try(npmle(nf,nz,plot="hist"),silent=TRUE),"try-error"))
