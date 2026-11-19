import os, tempfile
"""Create fixtures and direct-C expected results; never import the packages."""
import argparse, hashlib, json, sys, subprocess
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from raw_ctypes import Direct, POLICIES
from fixtures import simulate, save_arrays, encode, read_json, write_json
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--bunch", type=Path, required=True)
p.add_argument("--rd", type=Path)
p.add_argument("--root", type=Path, default=Path(os.environ.get("SPECPKG_BUILD", str(Path(tempfile.gettempdir())/"spectral-build"))) / 'xlang/fixtures')
p.add_argument("--entry", choices=["mass_gaussian","mass_voigt","mass_gaussian_jump","mass_joint_voigt","bmean","rd","npmle","postestimation"])
a = p.parse_args(); a.root.mkdir(parents=True, exist_ok=True)
if not a.entry:
    cases = []
    for entry in ("mass_gaussian","mass_voigt","mass_gaussian_jump","mass_joint_voigt","bmean","npmle","postestimation"):
        subprocess.run([sys.executable,"-I","-B",__file__,"--bunch",str(a.bunch),"--root",str(a.root),"--entry",entry],check=True)
        cases.extend(read_json(a.root/(entry+".json.gz"))["cases"])
    assert any(c["expected"]["status"] > 0 for c in cases), "missing refusal fixture"
    assert any(float.fromhex(c["expected"]["scalars"].get("tier", {"hex":"0x0p+0"})["hex"]) > 1 for c in cases), "missing fallback fixture"
    manifest={"format":"spectral-xlang-v1","byte_order":"little","float":"IEEE754-binary64","cases":cases,
              "libraries":{str(p):{"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in (a.bunch,)}}
    write_json(a.root/"manifest.json.gz",manifest)
    for entry in ("mass_gaussian","mass_voigt","mass_gaussian_jump","mass_joint_voigt","bmean","npmle","postestimation"):
        (a.root/(entry+".json.gz")).unlink()
    print("saved",len(cases),"cases",flush=True)
    raise SystemExit(0)
if a.entry == "postestimation":
    subprocess.run([sys.executable,"-I","-B",str(Path(__file__).resolve().parents[1]/"postestimation/generate.py"),"--library",str(a.bunch),"--root",str(a.root)],check=True)
    raise SystemExit(0)
if a.entry == "npmle":
    subprocess.run([sys.executable,"-I","-B",str(Path(__file__).resolve().parents[1]/"npmle/generate.py"),"--library",str(a.bunch),"--root",str(a.root)],check=True)
    (a.root/"manifest.json.gz").replace(a.root/"npmle.json.gz")
    raise SystemExit(0)
libs = {a.entry=="rd": Direct(a.rd if a.entry=="rd" else a.bunch, rd=a.entry=="rd")}
cases = []
for entry in (a.entry,):
    isrd = entry == "rd"
    policies = POLICIES if entry.startswith("mass") else [None]
    for shape in ("base", "skew", "refusal", "fallback"):
        arrays = simulate(1200, seed=17, shape=shape, rd=isrd)
        if not isrd and shape == "refusal":
            arrays["z"][:] = 0
        if not isrd and shape == "fallback":
            arrays["z"] += 0.1
        if isrd and shape == "refusal":
            arrays["d"][:] = 0
        if shape == "skew" and entry not in ("mass_gaussian_jump","mass_joint_voigt"):
            arrays["weights"] = np.ascontiguousarray(1.0 + (np.arange(1200)%3) * (1.0 if isrd else 0.25))
        cols = save_arrays(a.root, entry+"_"+shape, arrays)
        for policy in policies:
            options = {"seed": 1}
            if policy is not None: options["policy"] = policy
            case = {"id": entry+"_"+shape+"_"+str(policy), "entry": entry, "options": options, "columns": cols}
            status, raw = libs[isrd].run(entry, arrays, options)
            raw.pop("threads_used", None)
            case["expected"] = encode(status, raw)
            cases.append(case)
            print(case["id"], "status", status, "tier", raw.get("tier"), flush=True)
    # No large input files: each runner regenerates this deterministic recipe.
    recipe = {"n": 1000000, "seed": 23, "rd": isrd}
    arrays = simulate(**recipe)
    options = {"seed": 1}
    status, raw = libs[isrd].run(entry, arrays, options)
    raw.pop("threads_used", None)
    cases.append({"id": entry+"_large", "entry": entry, "options": options, "recipe": recipe, "expected": encode(status, raw)})
    print(entry, "large done", flush=True)
manifest = {"format": "spectral-xlang-v1", "byte_order": "little", "float": "IEEE754-binary64", "cases": cases,
            "libraries": {str(p): {"sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in (a.bunch,)},
            "build_info": [x.f("build_info")().decode() for x in libs.values()]}
write_json(a.root/(a.entry+".json.gz"),manifest)
print("saved", len(cases), "cases", flush=True)
