import os, tempfile
"""Install wheels in an isolated venv, compare all outputs and thread counts."""
import argparse, json, os, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
p = argparse.ArgumentParser()
p.add_argument("--root", type=Path, default=Path(os.environ.get("SPECPKG_BUILD", str(Path(tempfile.gettempdir())/"spectral-build"))))
p.add_argument("--worker", choices=["bunch","rd","both"])
p.add_argument("--tolerance", action="store_true", help="Cross-platform selected-answer contract; report all diagnostic differences")
a = p.parse_args()
if not a.worker:
    venv = a.root/"py/xlang-venv"
    env = dict(os.environ, UV_CACHE_DIR=str(a.root/"uvcache"), PYTHONDONTWRITEBYTECODE="1")
    subprocess.run(["uv", "venv", "--clear", str(venv)], env=env, check=True)
    python = venv/("Scripts/python.exe" if os.name == "nt" else "bin/python")
    wheels = sorted((a.root/"py/dist").glob("*.whl"))
    assert len(wheels) == 1, wheels
    subprocess.run(["uv", "pip", "install", "--python", str(python), *map(str,wheels), "pytest"], env=env, check=True)
    for family in ("both",):
        subprocess.run([str(python), "-I", "-B", __file__, "--root", str(a.root), "--worker",family]+(["--tolerance"] if a.tolerance else []), env=env, check=True)
    receipts=[json.loads((a.root/("py/xlang_"+f+".json")).read_text()) for f in ("both",)]
    total={k:sum(r[k] for r in receipts) for k in ("cases","calls","compared_values")}
    total.update(result="PASS",threads=[1,2,4,8,0],private_openmp=True)
    (a.root/"py/xlang_result.json").write_text(json.dumps(total,indent=2)+"\n")
    print(total)
else:
    from fixtures import load_arrays, public, encode, read_json
    root = a.root/"xlang/fixtures"
    manifest = read_json(root/"manifest.json.gz")
    manifest["cases"]=[c for c in manifest["cases"] if a.worker=="both" or (c["entry"]=="rd") == (a.worker=="rd")]
    from compare import CrossReport
    report = CrossReport(a.root/'cross-os/wheels') if a.tolerance else None
    values = 0
    bitwise = 0
    maxrel = 0.0
    for case in manifest["cases"]:
        arrays = load_arrays(root, case)
        expected = case["expected"]
        first_actual = None
        for threads in (1,2,4,8,0):
            result = public(case, arrays, threads)
            import specbunch
            default = specbunch.default_threads
            package = specbunch
            serial = "openmp=off" in package._core.fn("build_info")().decode()
            assert result.raw["threads_used"] == (1 if serial else (threads or default()))
            result.raw.pop("threads_used")
            actual = encode(result.status_code if case["entry"] in ("npmle","nonbuncher_mean","eti") else result.status, result.raw)
            from compare import compare
            stats = compare(actual, expected, exact=not a.tolerance, report=report, fixture=case['id'])
            if first_actual is not None:
                compare(actual, first_actual, exact=True)
            else:
                first_actual = actual
            bitwise += stats['bitwise']
            maxrel = max(maxrel, stats['maxrel'])
            values += 1 + len(expected["scalars"]) + len(expected["strings"]) + sum(map(len,expected["vectors"].values()))
        print(case["id"], "1/2/4/8/default PASS", flush=True)
    if report is not None: report.assert_pass()
    receipt = {"cases": len(manifest["cases"]), "calls": 5*len(manifest["cases"]), "compared_values": values,
               "threads": [1,2,4,8,0], "bitwise": bitwise, "maxrel": maxrel, "result": "PASS"}
    (a.root/("py/xlang_"+a.worker+".json")).write_text(json.dumps(receipt, indent=2)+"\n")
    print(receipt)
