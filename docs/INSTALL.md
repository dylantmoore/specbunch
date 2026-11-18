# Installation

Python: `python bunch/python/tools/build.py --root /path/to/build`, then install the wheel from `/path/to/build/py/dist`.

R: run `R CMD build /path/to/repo/bunch/r` from an external build directory, then `R CMD INSTALL` the resulting archive.

Julia: run `python tools/ci.py julia --root /path/to/build`.

Stata: run `python tools/ci.py stata --root /path/to/build --platform macarm64` (choose the platform matching Stata).
