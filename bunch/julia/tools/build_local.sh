#!/bin/sh
set -eu
R=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
B="${SB_JULIA_BUILD:-${SPECPKG_BUILD:-${TMPDIR:-/tmp}/spectral-build}/julia}"
export JULIA_DEPOT_PATH="$B/depot"
export TMPDIR="$B/tmp"
mkdir -p "$TMPDIR"
for family in bunch rd; do
    cmake -S "$R/$family/core" -B "$B/build/$family" -DCMAKE_BUILD_TYPE=Release
    cmake --build "$B/build/$family" --parallel 3
done
julia --startup-file=no "$R/bunch/julia/tools/install_local.jl"
