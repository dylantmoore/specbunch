args buildroot
if `"`buildroot'"' == "" local buildroot : environment SPECPKG_BUILD
clear all
set more off
sysdir set PLUS "`buildroot'/plus"
sysdir set PERSONAL "`buildroot'/personal"
net install specbunch, from("`buildroot'/net/specbunch") replace
set obs 600
set seed 12345
generate double z=round(rnormal(),1/1024)
generate double x=2*(z>0)+rnormal()*.1
generate double d=z>0
generate double wt=1
replace z=.a in 5
replace x=. in 6
replace wt=0 in 7
replace wt=. in 8

specbunch buncher z x [pw=wt] if _n>2 in 1/550, threads(1)
assert e(n)==544
local status=e(core_status)
scalar expected_mu=e(mu)
scalar expected_se=e(se)
preserve
keep in 1/550
drop in 1/2
drop if missing(z,x,wt) | wt<=0
specbunch buncher z x [pw=wt], threads(1)
assert e(core_status)==`status'
assert e(mu)==expected_mu
assert e(se)==expected_se
restore
generate double shifted=z+8
specbunch mass shifted, kink(8) threads(1)
local bits "`e(scalar_bits)'"
specbunch mass z, threads(1)
assert "`e(scalar_bits)'"=="`bits'"
replace shifted=z-8
specbunch mass shifted, kink(-8) threads(1)
assert "`e(scalar_bits)'"=="`bits'"
assert "`e(policy)'"=="tiered"
assert colsof(e(windows)) > 10
capture specbunch mass z [pw=wt], method(joint) friction(gaussian) threads(1)
assert _rc==498
specbunch mass z in 1/10, threads(1)
assert "`e(str_reason)'"=="insufficient_data"
capture specbunch mass z, seed(-1) threads(1)
assert _rc==498
capture specbunch mass z, policy(bad)
assert _rc==198
specbunch buncher z x, threads(1) seed(18446744073709551615) draws(19)
assert "`e(cmd)'"=="specbunch"
foreach t in 1 3 16 {
    specbunch buncher z x, draws(19) threads(`t')
    assert e(threads_used) == `t'
}
specbunch buncher z x, draws(19)
local expected : environment SPECTRAL_EXPECT_DEFAULT
if "`expected'" == "" local expected 4
assert e(threads_used) == real("`expected'")
clear
set obs 1500
set seed 37
generate double z=cond(runiform()<.2,0,rnormal())+.04*rnormal()
generate double weights=1+mod(_n,3)*.25
specbunch mass z, method(joint) friction(voigt) threads(1)
assert e(core_status)==0
assert "`e(method)'"=="joint" & "`e(friction)'"=="voigt"
assert rowsof(e(transform))==401 & colsof(e(transform))==8
assert e(a)>=0 & e(s2)>=0 & e(sigma_eq)>0
specbunch mass z [pw=weights], threads(1)
assert e(core_status)==0
assert rowsof(e(transform))==401
assert colsof(e(transform))==5
assert colnumb(e(transform),"cosine_fit")<.
assert e(n_points)>=11
specbunch npmle, threads(1) gridpoints(101) lambdapoints(5) refinepoints(3)
assert "`e(npmle_status)'"=="converged"
assert e(npmle_kkt)<.0000001
assert colsof(e(npmle_curves))==6
assert colsof(e(npmle_histogram))==7
assert e(npmle_B)==e(B)
quietly count if e(sample)
local in_sample=r(N)
graph drop _all
specbunch npmle, threads(1) gridpoints(101) lambdapoints(3) refinepoints(3) name(sbnpmle, replace)
quietly graph dir
assert strpos(" `r(list)' "," sbnpmle ")>0
quietly count if e(sample)
assert r(N)==`in_sample'
assert "`e(npmle_status)'"=="converged"
specbunch npmle, threads(1) gridpoints(101) lambdapoints(3) refinepoints(3) plot(cdf) name(sbcdf)
quietly graph dir
assert strpos(" `r(list)' "," sbcdf ")>0
graph drop _all
specbunch npmle, threads(1) gridpoints(101) lambdapoints(3) refinepoints(3) noplot
quietly graph dir
assert strtrim("`r(list)'")==""
capture specbunch npmle, plot(hist)
assert _rc==198
specbunch mass z in 1/5, threads(1)
capture specbunch npmle
assert _rc==498
display "UNIT_SPECBUNCH_PASS"
exit, clear
