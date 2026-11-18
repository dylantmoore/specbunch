local buildroot : environment SPECPKG_BUILD
clear all
set more off
sysdir set PLUS "`buildroot'/stata/plus"
sysdir set PERSONAL "`buildroot'/stata/personal"
net install specbunch, from("`buildroot'/stata/net/specbunch") replace
net install specrd, from("`buildroot'/stata/net/specrd") replace
set obs 600
set seed 45
gen double z=rnormal()
gen double d=z>0
gen double y=2*d+.1*rnormal()
foreach threads in 1 2 4 8 {
    specrd y d z, threads(`threads')
    assert e(n)==600
    specbunch buncher z y, threads(`threads')
    assert e(n)==600
}
display "COEXIST_REVERSE_PASS"
exit, clear
