# Cross-language consistency QA

All inputs are simulated. `generate.py --bunch LIB --root FIXTURES`
creates the spectral-xlang-v1 manifest and small raw column files. Expected
results come directly from the public C ABI. Large samples are stored as recipes,
not data. The generator uses ordered dyadic arithmetic and counter-based uint64
hashing, so input bytes are portable across platforms.

Each number retains its IEEE754 bits and hexadecimal value. Compare every field,
including status, decisions, flags, tiers, strings, vectors and key sets.
`run_python.py --root BUILD` installs wheels and tests 1/2/4/8/default threads.
Its default is bitwise comparison; `--tolerance` allows relative error <= 1e-10
in numerical fields while requiring exact discrete outputs. Bitwise counts are
reported in both modes. Thread counts are checked separately from fixture outputs.

`run_r.R`, `run_julia.jl` and `run_stata.py` exercise each language's public API.
The Stata runner requires `STATA_EXECUTABLE`; the CI `tests/no_stata` harness
instead loads the same compiled plugin and supplies SPI callbacks in memory.
Fixtures are generated once on Linux in CI and downloaded by platform jobs.
macOS arm64 additionally compares its wheels to a native C-ABI oracle exactly;
Linux-to-macOS comparisons use tolerance because the system math libraries differ.

Deployment audits are platform-specific. `tests/test_threads.cpp`
checks detection, environment overrides and restoration. Outputs belong outside
the checkout. Fixtures and CI artifacts have seven-day retention; no collectors
or recurring backup processes are installed.
