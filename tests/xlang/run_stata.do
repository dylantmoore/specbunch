* Public Stata cross-language runner. Generate the exact fixture plan first:
* python -B run_stata.py --prepare-only (requires numpy).
args buildroot
if `"`buildroot'"' == "" local buildroot : environment SPECPKG_BUILD
do "`buildroot'/xlang/plan.do"
