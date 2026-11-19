# ETI shoulder validation

All inputs are simulated. `fixture.cpp` checks clamped linear basis normalization,
checks Gaussian convolution against independent Simpson integration, tests
zero/negative-shoulder statuses, and emits the two-basis C ABI fixture.
`expected.csv` is new; no existing fixture is replaced. Its copies in the
Python/R/Julia test directories are byte-identical. Wrapper checks use the
package's cross-platform tolerance, while the existing fixture sweeps retain
bitwise checks within a build.

Build and run `sb_test_eti_shoulders` from the core CMake project. Python's
`test_eti_shoulders.py`, R's `eti_shoulders.R`, Julia's final testset, Stata's
`test_eti_shoulders.do NET_DIRECTORY FIXTURE_DIRECTORY`, and
`no_stata.py HOST_LIBRARY PLUGIN expected.csv` consume these fixtures.

`mc.cpp` is a standalone C++20 driver linked to libspecbunch and an OpenMP runtime.
It accepts `REPS THREADS OUTPUT_CSV CELLS [FIRST_CELL]`; the complete experiment is
`100 8 output.csv 32`. It fixes first-stage mass/friction and the transform band
at known values to isolate shoulder error. Each cell is sample size × friction
width × jump × friction law; all methods use the same sample. There are no CIs.
Numerical settings are explicit in the source. Standard-library RNG results are
deterministic within a platform/build; estimates are invariant to worker count.
The final CSV is sorted by cell and repetition regardless of scheduling.

Example on Linux, with an existing core build at `$BUILD/core` (put `$BUILD`
under `~/scratch-local/spectral-bunching/pkg_build`):

```sh
c++ -std=c++20 -O3 -ffp-contract=off -fopenmp \
  -Ibunch/core/include -Icommon tests/eti_shoulders/mc.cpp \
  -L"$BUILD/core" -lspecbunch -Wl,-rpath,"$BUILD/core" -o "$BUILD/eti_mc"
"$BUILD/eti_mc" 100 8 "$BUILD/mc.csv" 32 > "$BUILD/mc.log" 2>&1
python3 tests/eti_shoulders/report.py "$BUILD/mc.csv" "$BUILD/mc.log" "$BUILD/ETI_SHOULDERS_MC.md"
```

The driver retains only one sample per worker, one fit at a time, and compact
scalar records. It creates no raw-data files, database snapshots, or ongoing
collector. CSV totals and cell counts are checked before generating the report.
