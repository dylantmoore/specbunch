# specbunch

Bunching estimators for a kink with optimization frictions, in Python, R, Julia and Stata. All four call the same deterministic C++ core.

`specbunch` estimates the bunching mass B at a kink when bunchers are spread out by a friction (reporting or adjustment error, say), the mean characteristic of bunchers and of nonbunchers at the kink, the elasticity of taxable income (ETI), and a counterfactual density by NPMLE. The key idea: friction convolves the latent distribution, which multiplies its characteristic function by the friction's characteristic function, so an atom of mass B at the kink appears as `B * phi(t)` in the empirical cosine transform at frequencies where the smooth part of the density has died out; B and the friction scale are fitted there directly.

- **Using it with an AI assistant?** Point it to [SKILL.md](SKILL.md): a compact guide to the API, the output fields and flags, and the recommended workflow.
- **Every option, field, flag and refusal reason:** [docs/REFERENCE.md](docs/REFERENCE.md).
- **Platforms and reproducibility:** [docs/PLATFORMS.md](docs/PLATFORMS.md).
- **Installation details and troubleshooting:** [docs/INSTALL.md](docs/INSTALL.md).

## Status and testing

This is pre-release research software. Please read this before relying on it.

- **Tested** on macOS Apple Silicon in all four languages, including real Stata (StataNow), with extensive cross-language tests that require bitwise-identical results from Python, R, Julia and Stata on the same machine.
- **Continuous integration** on GitHub-hosted runners builds and tests the Python package on Linux x86_64 and aarch64, macOS Intel and Apple Silicon, and Windows x86_64, and the R and Julia packages on Linux, macOS and Windows.
- **Stata on Windows, Linux and Intel Mac:** the plugins are built and tested in CI only through a harness that imitates Stata's plugin interface, not in real Stata.
- **Not yet tested on users' own machines.** Windows is overall the least tested platform.
- **Reproducibility:** on one platform and build, results are bitwise identical at any thread count. Across operating systems the selected answer (status, tier, flags, window, B, sigma, means) agreed on all test fixtures, but diagnostics such as fit p-values can differ in the last digits, and a sample on a knife edge between two windows could select a different one. Report the platform with published results. See [docs/PLATFORMS.md](docs/PLATFORMS.md).
- **Standard errors are model-based** and conditional on the selected window. They are not confidence intervals.
- **Estimators and defaults may change** before the accompanying paper is final.
- **Please report problems**, installation failures included, at <https://github.com/dylantmoore/specbunch/issues>, with your operating system, language version and the output of `version()` (Stata: `which specbunch`).

## Installation

Nothing is on PyPI, CRAN, the Julia General registry or SSC yet. Every route below builds from source, so you need a C++20 compiler and an internet connection during the build (it downloads two pinned, checksummed dependencies: the HiGHS LP solver and, on macOS and Windows, the LLVM OpenMP runtime; the R package bundles its sources and downloads nothing). Prebuilt binaries (Python wheels, Stata plugins) are planned, so that no compiler will be needed.

### Prerequisites

| OS | What to install |
|---|---|
| macOS | Xcode Command Line Tools: `xcode-select --install`. For Julia and Stata also CMake ≥ 3.20 (`brew install cmake` or <https://cmake.org/download/>) and Python 3 (to run the Stata build script). |
| Linux | GCC or Clang with C++20 support, `git`, CMake ≥ 3.20 (Debian/Ubuntu: `sudo apt install build-essential git cmake`); R also needs `r-base-dev`. Python wheels and the Stata plugin use the system `libgomp`. |
| Windows | R: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version. Python, Julia and Stata: Visual Studio Build Tools (C++ workload), LLVM (`clang-cl`), CMake, Ninja and Git; see [docs/INSTALL.md](docs/INSTALL.md#windows). |

The Python build fetches CMake and Ninja itself if they are missing.

### Python (3.10+)

```bash
python -m pip install "specbunch[plot] @ git+https://github.com/dylantmoore/specbunch#subdirectory=bunch/python"
```

The `[plot]` extra installs matplotlib for the default plots; for a minimal install (NumPy only) use `"git+https://github.com/dylantmoore/specbunch#subdirectory=bunch/python"`. Compiling takes a few minutes.

### R (4.2+)

```r
install.packages("remotes")
remotes::install_github("dylantmoore/specbunch", subdir = "bunch/r")
```

macOS needs the Xcode Command Line Tools and Windows needs Rtools (see above). Compiling takes several minutes. With CRAN's macOS build of R, which has no OpenMP, the package runs single-threaded; results are identical, only slower.

### Julia (1.9+)

The binary package `SpecBunch_jll` is not registered yet, so for now you build the core library with CMake and point a small local stand-in package at it. In a shell, from a directory of your choice:

```bash
git clone https://github.com/dylantmoore/specbunch
cmake -S specbunch/bunch/core -B specbunch-build -DCMAKE_BUILD_TYPE=Release
cmake --build specbunch-build --config Release --parallel
```

Then start `julia` in the same directory (activate the environment you want it in, or stay in the default one) and run:

```julia
using Pkg
repo  = abspath("specbunch")           # the git clone
build = abspath("specbunch-build")     # the CMake build directory
lib = joinpath(build, Sys.iswindows() ? "specbunch.dll" :
                      Sys.isapple()   ? "libspecbunch.dylib" : "libspecbunch.so")
isfile(lib) || error("core library not found: $lib")
# Local stand-in for the not-yet-registered SpecBunch_jll binary package.
jll = joinpath(build, "SpecBunch_jll"); mkpath(joinpath(jll, "src"))
write(joinpath(jll, "Project.toml"),
      "name = \"SpecBunch_jll\"\nuuid = \"7035fe15-a909-5e98-9b00-51961afcec2a\"\nversion = \"0.1.0\"\n")
write(joinpath(jll, "src", "SpecBunch_jll.jl"),
      "module SpecBunch_jll\nexport libspecbunch\nconst libspecbunch = $(repr(lib))\nend\n")
Pkg.develop([PackageSpec(path=jll), PackageSpec(path=joinpath(repo, "bunch", "julia", "SpecBunch"))])
```

Keep both directories: the package loads the library from `specbunch-build`. Plotting is a package extension that `using Plots` activates.

### Stata (64-bit)

Build the plugin for your platform (`macarm64` for Apple Silicon, `macintel64` for Intel Mac, `linux64` or `win64`), then install from the local folder. In a shell:

```bash
git clone https://github.com/dylantmoore/specbunch
python3 specbunch/bunch/stata/build.py --platform macarm64 --build-root "$PWD/specbunch-stata"
```

The `--build-root` path must be absolute (`$PWD/...` does that; on Windows use `"%CD%\specbunch-stata"` in Command Prompt or `"$PWD\specbunch-stata"` in PowerShell). Then in Stata, with the full path to that folder:

```stata
net install specbunch, from("/full/path/to/specbunch-stata/net/specbunch") replace
help specbunch
```

macOS Gatekeeper: a plugin you build yourself is not quarantined and loads normally. A `.plugin` file downloaded through a browser is unsigned, and macOS may refuse to load it; clear the quarantine flag with `xattr -d com.apple.quarantine specbunch_macarm64.plugin` in the folder where it was installed (`sysdir` shows the PLUS directory; the file is under `s/`).

### Check the installation

| Language | Command |
|---|---|
| Python | `python -c "import specbunch; print(specbunch.version())"` |
| R | `library(specbunch); version()` |
| Julia | `using SpecBunch; SpecBunch.version()` |
| Stata | `which specbunch` |

**Verified** on 2026-11-19 on macOS Apple Silicon only: each route above was run literally against this public repository in a clean environment (fresh Python venv, fresh R library, fresh Julia depot, isolated Stata PLUS directory with StataNow in batch mode), followed by a mass estimate. The Linux and Windows instructions follow the CI configuration but have not been run by hand.

## Quick start (Python)

```python
import numpy as np, specbunch

z = np.log(income) - np.log(kink)             # log distance to the kink; full sample, no trimming
fit = specbunch.mass(z)                       # Gaussian friction, cosine method
print(fit.B, fit.sigma, fit.tier, fit.flags)  # estimate, selection tier, diagnostic flags
fit.windows                                   # every candidate window (a DataFrame if pandas is installed)
voigt = specbunch.mass(z, friction="voigt")   # report Gaussian and Voigt as a pair
```

R, Julia and Stata use the same names; [SKILL.md](SKILL.md) has a per-language table.

## Methods

Methods: `cosine` (default), `joint`. Friction: `gaussian` (default), `voigt`. Window rules: `tiered` (default), `first`. Every fit returns ranked sensitivity tables; coherent estimates carry diagnostic flags.

The mass estimator is weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform's own variance, solved by Gauss–Newton/Newton. The joint method also fits sine; the companions are log-linear OLS/LAD.

`npmle` adds fixed-mass, fixed-friction spline postestimation, with a BIC-selected penalty, histogram shares and latent/observed/empirical curves. Each language wrapper draws default plots: `mass` shows the empirical transform with its fitted curve and selected window, and `npmle` shows the fitted distributions; pass `plot=False` (Stata `noplot`) to skip them.

Postestimation also provides `nonbuncher_mean` (a local characteristic mean with the buncher mean fixed) and ETI by the Saez trapezoid with continuous completion. See [docs/REFERENCE.md](docs/REFERENCE.md) for language-specific syntax, conditional SE interpretation and flags, and the `VENDOR.md` files for source provenance. Build artifacts must be placed outside the source tree.

## Citation

If you use `specbunch`, please cite:

> Moore, Dylan (2026). "Spectral Deconvolution: Applications to Bunching and Regression Discontinuity Designs." Working paper.

```bibtex
@unpublished{moore2026spectral,
  author = {Moore, Dylan},
  title  = {Spectral Deconvolution: Applications to Bunching and Regression Discontinuity Designs},
  note   = {Working paper},
  year   = {2026}
}
```
