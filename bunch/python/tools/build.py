import os, tempfile
"""Stage self-contained sources outside Drive, then build sdist + wheel with uv.
Usage: python -B tools/build.py [--stage-only] [--root /external/build/root]
"""
import argparse, os, platform, re, shutil, subprocess, sys
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument("--stage-only", action="store_true")
p.add_argument("--arch", choices=["arm64", "x86_64"], help="macOS target architecture (default: ARCHFLAGS or host)")
p.add_argument("--clean", action="store_true", help="Remove this package's external build tree before rebuilding")
p.add_argument("--root", type=Path, default=Path(os.environ.get("SPECPKG_BUILD", str(Path(tempfile.gettempdir())/"spectral-build"))))
a = p.parse_args()
src = Path(__file__).resolve().parents[1]
part = src.parent.name
name = "specbunch" if part == "bunch" else "specrd"
root = a.root.expanduser().resolve()
if "CloudStorage" in str(root) or "GoogleDrive" in str(root):
    raise SystemExit("build root must be outside Drive")
stage = root/"py/stage"/name
if stage.exists():
    shutil.rmtree(stage)
shutil.copytree(src, stage, ignore=shutil.ignore_patterns("__pycache__", ".pytest_cache", "native", "dist", "build"))
packages = src.parents[1]
# Read legal files explicitly: macOS File Provider fast-copy can copy an
# unhydrated placeholder as an empty file. Never ship empty licence texts.
for rel in [Path(x) for x in ("LICENSE", "NOTICE", "THIRD_PARTY_NOTICES.md", "CITATION.cff")] + [Path("licenses")/x.name for x in (packages/"licenses").iterdir()]:
    data=(packages/rel).read_bytes()
    if not data: raise SystemExit("Empty legal file: "+str(rel))
    (stage/rel).parent.mkdir(parents=True,exist_ok=True)
    (stage/rel).write_bytes(data)

shutil.copytree(packages/part/"core", stage/"native"/part/"core")
shutil.copytree(packages/"common", stage/"native/common")
subprocess.run([sys.executable, "-B", str(packages/"tools/source_archive.py"),
                "--repo", str(packages), "--output", str(stage/"src/specbunch/licenses/corresponding-source.tar.gz"),
                "--cache", str(root/"source-cache")], check=True)
print(stage, flush=True)
if not a.stage_only:
    env = dict(os.environ, UV_CACHE_DIR=str(root/"uvcache"), PYTHONDONTWRITEBYTECODE="1", CMAKE_BUILD_PARALLEL_LEVEL="4")
    tmp = root/"py/tmp"; tmp.mkdir(parents=True, exist_ok=True)
    env["TMPDIR"] = str(tmp)
    settings=[]
    if sys.platform == "darwin":
        requested=re.findall(r"-arch\s+(\S+)",env.get("ARCHFLAGS",""))
        if len(set(requested))>1:
            raise SystemExit("Build arm64 and x86_64 separately to preserve their deployment floors")
        arch=a.arch or (requested[0] if requested else platform.machine())
        if arch not in ("arm64","x86_64"):
            raise SystemExit("Unsupported macOS architecture: "+arch)
        target="11.0" if arch=="arm64" else "10.15"
        env["ARCHFLAGS"]="-arch "+arch
        env["MACOSX_DEPLOYMENT_TARGET"]=target
        settings=["-C", "cmake.define.CMAKE_OSX_ARCHITECTURES="+arch,
                  "-C", "cmake.define.CMAKE_OSX_DEPLOYMENT_TARGET="+target]
    build=root/"py/build"/name
    if a.clean and build.exists():
        shutil.rmtree(build)
    subprocess.run(["uv", "build", str(stage), "--out-dir", str(root/"py/dist"),
                    "-C", "build-dir="+str(build), *settings], env=env, check=True)
