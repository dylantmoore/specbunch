
## Spline postestimation
The authored npmle implementation adapts polynomial B-spline convolution and active-set simplex Newton ideas from spectral-bunching/cpp/part1/src/spl/spl.hpp and eti/fcf.hpp. The fixed-parameter likelihood, second-difference penalty and BIC trace criterion follow nta_slides/talk/prep_npmle.py and cjms-estimators/src/spline_deconvolution.jl. Private input loading and plotting code are not used.

## Postestimation provenance

Non-vendored `src/postestimation.cpp` implements the fixed-buncher-mean triangular
local polynomial described in `~/scratch-local/spectral-bunching/cefmix/`
(`METHOD.md`, `NOTES.md`, `cefmix.hpp`, `capi.cpp`, `run_lp.py`, `run_fix.py`).
The source is read-only reference material; empirical input files were not used.
ETI implements the CF-native midpoint/jump and Saez/C0 closed forms in
`spectral-bunching/cpp/part1/src/eti/eti.hpp` and `etilib.cpp`, with the scope
in `paper/outline_2026-11-11.md`, Appendix E. `bfc.hpp` was reviewed but its
counterfactual optimisation is not included. No forward-CF or C1/C2 rules are
exposed. Public names are `saez` and `continuous`. Existing vendored files remain
byte-identical. Validation compiles the reference headers separately on simulated data.
