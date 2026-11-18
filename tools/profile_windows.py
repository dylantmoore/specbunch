"""Instrument a disposable core copy; never alter the shipped source tree."""
import argparse, json, os, shutil, subprocess, sys, time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--worker",nargs=3);p.add_argument("--library");a=p.parse_args();repo=Path(__file__).resolve().parents[1]
if a.worker:
    sys.path.insert(0,str(repo/"tests/xlang"))
    import numpy as np
    from raw_ctypes import Direct
    from fixtures import simulate
    entry,kind,mode=a.worker
    if kind=="slow_test":
        rng=np.random.default_rng(42);z=np.where(rng.random(1000)<0.3,0,rng.normal(size=1000))+rng.normal(0,0.04,1000);arrays={"z":z}
    else:arrays=simulate(1200,seed=17)
    lib=Direct(a.library);start=time.perf_counter();status,raw=lib.run(entry,arrays,{"seed":1,"check_all_windows":int(mode)},threads=1)
    print(json.dumps(dict(entry=entry,kind=kind,check_all_windows=int(mode),seconds=time.perf_counter()-start,status=status,scalars={k:v for k,v in raw.items() if isinstance(v,(float,int))})),flush=True)
    raise SystemExit()
root=a.root.resolve();src=root/"source";evidence=root/"evidence";evidence.mkdir(parents=True,exist_ok=True)
for part in ("bunch/core","common"):
    shutil.copytree(repo/part,src/part,dirs_exist_ok=True)
f=src/"bunch/core/src/distribution_screen.hpp";s=f.read_text();s=s.replace("namespace sbdistribution {",r'''namespace sbdistribution {
struct AuditTimer {const char* name; double B, sigma; std::chrono::steady_clock::time_point start;
 AuditTimer(const char* n,double b=0,double sig=0):name(n),B(b),sigma(sig),start(std::chrono::steady_clock::now()){}
 ~AuditTimer(){double t=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::fprintf(stderr,"PROFILE %s %.9f B=%.17g sigma=%.17g\n",name,t,B,sigma);std::fflush(stderr);}
};
''')
s="#include <chrono>\n#include <cstdio>\n"+s
for start,label,extra in [("inline void herm_eig(","herm_eig",""),("inline LPRes linprog_highs(","linprog_highs",""),("inline double lower_polynomial(","lower_polynomial",""),("inline std::string compat_screen(","compat_screen",", B, sigma"),("  std::string check(double B, double sigma)","check",", B, sigma"),("inline int check_all(","check_all",", B, sigma")]:
    pos=s.index("{",s.index(start))+1;s=s[:pos]+f'\n AuditTimer audit_timer("{label}"{extra});\n'+s[pos:]
f.write_text(s)
env=os.environ.copy();env.update(CC="clang-cl",CXX="clang-cl",CMAKE_GENERATOR="Ninja")
build=root/"build"
subprocess.run(["cmake","-S",str(src/"bunch/core"),"-B",str(build),"-DCMAKE_BUILD_TYPE=Release"],env=env,check=True)
subprocess.run(["cmake","--build",str(build),"--parallel","3"],env=env,check=True)
lib=build/"specbunch.dll"
for entry,kind,mode,limit in [("mass_gaussian","slow_test",0,180),("mass_gaussian","slow_test",2,420)]+[(e,"base",1,180) for e in ("mass_gaussian","mass_voigt","mass_gaussian_jump")]:
    name=f"{entry}-{kind}-{mode}";cmd=[sys.executable,"-B",__file__,"--root",str(root),"--library",str(lib),"--worker",entry,kind,str(mode)]
    with (evidence/(name+".out")).open("w") as out,(evidence/(name+".trace")).open("w") as err:
        try:
            r=subprocess.run(cmd,stdout=out,stderr=err,timeout=limit);result=dict(case=name,returncode=r.returncode)
        except subprocess.TimeoutExpired:result=dict(case=name,timeout=limit)
    with (evidence/"profile-status.jsonl").open("a") as f:print(json.dumps(result),file=f)
    print(result,flush=True)
