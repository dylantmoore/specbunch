# Run in a VS x64 developer shell with clang-cl, CMake and Python3.
# Private static LLVM OpenMP uses the core pinned Windows patch.
$env:CC="clang-cl"
$env:CXX="clang-cl"
$env:CMAKE_GENERATOR="Ninja"
python -B "$PSScriptRoot/build.py" --platform win64
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
