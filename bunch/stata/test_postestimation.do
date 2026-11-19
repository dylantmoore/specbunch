version 14.1
set more off
args buildroot
adopath ++ "`buildroot'/specbunch/macarm64"
clear
set obs 5000
generate double z=-.5+(_n-1)/4999
generate double x=.4
program define fixed_mass, eclass
    tempname b
    matrix `b'= (.08)
    ereturn post `b'
    ereturn local cmd specbunch
    ereturn local method cosine
    ereturn local friction gaussian
    ereturn scalar core_status=0
    ereturn scalar B=.08
    ereturn scalar sigma=.03
    ereturn scalar a=0
    ereturn scalar kink=0
    ereturn scalar window_lo=12
    ereturn scalar window_hi=45
end
fixed_mass
estimates store mass_fit
specbunch nonbuncher_mean z x, massfit(mass_fit) mub(.4)
assert e(core_status)==0
assert abs(e(m_NB)-.4)<1e-10
assert rowsof(e(local_fit))==40
assert colnumb(e(local_fit),"fitted_mixture")<.
program define fixed_mean, eclass
    tempname b
    matrix `b'=(.4)
    ereturn post `b'
    ereturn local cmd specbunch
    ereturn local method buncher_mean
    ereturn scalar core_status=0
    ereturn scalar kink=0
    ereturn scalar mu=.4
    ereturn scalar se=.01
end
fixed_mean
estimates store mean_fit
specbunch nonbuncher_mean z x, massfit(mass_fit) buncherfit(mean_fit)
assert e(mu_b)==.4
assert e(se_mu_b)==.01
assert e(se)>0
specbunch nonbuncher_mean z x, massfit(mass_fit) mub(.4) degree(2)
assert abs(e(m_NB)-.4)<1e-10
specbunch eti z, massfit(mass_fit) netbelow(1.34) netabove(1)
assert e(core_status)==0
assert abs(e(saez_e)-e(saez_d)/ln(1.34))<1e-14
capture specbunch eti z, massfit(mass_fit) netbelow(1) netabove(1)
assert _rc!=0
display "POSTESTIMATION_STATA_PASS"
