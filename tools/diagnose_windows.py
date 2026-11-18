"""Read-only performance probes; each timed call uses a fresh process."""
import argparse, hashlib, json, os, platform, subprocess, sys, time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests/xlang"))
from fixtures import simulate, encode
from raw_ctypes import Direct
p=argparse.ArgumentParser();p.add_argument("--library",required=True);p.add_argument("--output",required=True);p.add_argument("--worker",nargs=3);a=p.parse_args()
if a.worker:
    entry,n,threads=a.worker
    lib=Direct(a.library);arrays=simulate(int(n),seed=17 if int(n)==1200 else 23)
    start=time.perf_counter(); status,raw=lib.run(entry,arrays,{"seed":1},threads=int(threads)); elapsed=time.perf_counter()-start
    used=raw.pop("threads_used",None)
    print(json.dumps(dict(entry=entry,n=int(n),threads=int(threads),threads_used=used,seconds=elapsed,status=status,tier=raw.get("tier"),build_info=lib.f("build_info")().decode(),sha256=hashlib.sha256(json.dumps(encode(status,raw),sort_keys=True).encode()).hexdigest())),flush=True)
else:
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as f:
        meta=dict(platform=platform.platform(),cpu_count=os.cpu_count(),processor=platform.processor())
        print(json.dumps(meta),file=f,flush=True);print(meta,flush=True)
        for mode in ("default","passive"):
            env=os.environ.copy()
            if mode=="passive":env.update(OMP_WAIT_POLICY="PASSIVE",KMP_BLOCKTIME="0")
            for entry,n in [(e,1200) for e in ("mass_gaussian","mass_voigt","mass_gaussian_jump")]+[("mass_gaussian",1000000),("mass_gaussian_jump",1000000)]:
                for threads in (1,2,4,8,0):
                    cmd=[sys.executable,"-B",__file__,"--library",a.library,"--output",a.output,"--worker",entry,str(n),str(threads)]
                    try:
                        result=subprocess.run(cmd,env=env,capture_output=True,text=True,timeout=100)
                        row=json.loads(result.stdout.strip()) if result.returncode==0 else dict(error=result.stderr[-3000:],returncode=result.returncode,entry=entry,n=n,threads=threads)
                    except subprocess.TimeoutExpired:row=dict(timeout=100,entry=entry,n=n,threads=threads)
                    row["wait_mode"]=mode
                    print(json.dumps(row),file=f,flush=True);print(json.dumps(row),flush=True)
