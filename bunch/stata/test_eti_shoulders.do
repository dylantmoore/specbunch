clear all
set more off
args net fixtures
adopath ++ "`net'"
program _specbunch_plugin, plugin using("`net'/specbunch_macarm64.plugin")
import delimited using "`fixtures'/expected.csv", clear varnames(1) asdouble
scalar expected_cubic=saez[1]
scalar expected_linear=saez[2]
clear
set obs 400
gen double z=(_n-1-200)/256
program fixed, eclass
 tempname b
 matrix `b'=(.08)
 ereturn post `b'
 ereturn local cmd specbunch
 ereturn local method cosine
 ereturn local friction gaussian
 ereturn scalar core_status=0
 ereturn scalar B=.08
 ereturn scalar sigma=.03
 ereturn scalar kink=0
 ereturn scalar window_lo=12
 ereturn scalar window_hi=45
 ereturn local running_variable z
 tempvar t
 gen byte `t'=1
 ereturn repost, esample(`t')
end
foreach basis in linear cubic {
 fixed
 specbunch npmle, basis(`basis') noplot threads(1) coreintervals(4) gridpoints(25) lambdapoints(3) refinepoints(3)
 scalar left=e(npmle_latent_left)
 scalar right=e(npmle_latent_right)
 specbunch eti z, netbelow(1.34) netabove(1) threads(1)
 assert abs(e(npmle_saez_e)-2*.08/(left+right)/log(1.34))<1e-14
 assert e(npmle_shoulder_left)==left
 assert abs(e(npmle_saez_e)-expected_`basis')<1e-8
 fixed
 specbunch npmle, basis(`basis') noplot threads(1) coreintervals(4) gridpoints(25) lambdapoints(3) refinepoints(3) netbelow(1.34) netabove(1)
 assert e(npmle_saez_e)<.
}
display "ETI_SHOULDERS_STATA_PASS"
exit, clear
