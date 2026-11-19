library(specbunch)
z <- ((0:399)-200)/256
fixture <- read.csv("eti_shoulders.csv")
f <- structure(list(B=.08,sigma=.03,a=0,friction="gaussian",kink=0,window_lo=12,window_hi=45,refused=FALSE),class="specbunch_mass")
for(b in c("linear","cubic")) {
 n <- npmle(f,z,basis=b,plot=FALSE,core_intervals=4,grid_points=25,lambda_points=3,refinement_points=3,threads=1,net_below=1.34,net_above=1)
 row <- fixture[fixture$basis==b,]
 stopifnot(abs(n$latent_left-row$left)<1e-8,abs(n$eti$npmle$saez$e-row$saez)<1e-8)
 stopifnot(n$eti$npmle$basis==b,n$eti$npmle$shoulder_left==n$latent_left,n$latent_left==n$coefficients[1])
 stopifnot(abs(n$eti$npmle$saez$e-2*f$B/(n$latent_left+n$latent_right)/log(1.34))<1e-14)
 f$npmle <- n
 stopifnot(eti(f,z,net_below=1.34,net_above=1)$npmle$saez$e==n$eti$npmle$saez$e)
}
