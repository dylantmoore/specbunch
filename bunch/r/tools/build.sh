#!/bin/sh
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
B=${SPECTRAL_R_BUILD:-"${SPECPKG_BUILD:-${TMPDIR:-/tmp}/spectral-build}/r"}
mkdir -p "$B/lib" "$B/tmp" "$B/logs"
export R_LIBS_USER="$B/lib" TMPDIR="$B/tmp"
export _R_CHECK_FORCE_SUGGESTS_=false _R_CHECK_CRAN_INCOMING_REMOTE_=false _R_CHECK_SYSTEM_CLOCK_=false
PATH="$B/bin:$PATH"; export PATH
cd "$B"
R CMD build "$here" > "logs/build_specbunch.log" 2>&1
R CMD INSTALL --library="$B/lib" specbunch_0.1.0.tar.gz > "logs/install_specbunch.log" 2>&1
R CMD check --as-cran specbunch_0.1.0.tar.gz > "logs/check_specbunch.log" 2>&1
