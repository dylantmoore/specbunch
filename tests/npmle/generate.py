"""Compact deterministic NPMLE oracle, independently through the C ABI."""
import argparse,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"xlang"))
from fixtures import simulate,save_arrays,encode,write_json
from raw_ctypes import Direct
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--root',type=Path,required=True);args=p.parse_args();args.root.mkdir(parents=True,exist_ok=True)
core=Direct(args.library);cases=[]
for name,n,method,weighted,kink in [('exact',600,'gaussian',False,0.),('weighted',2200,'gaussian',True,.125),('voigt',5000,'voigt',False,0.),('jump',50000,'gaussian_jump',False,0.)]:
 a=simulate(n,seed=29);a={'z':a['z']+kink}
 if weighted:a['weights']=1+(np.arange(n)%3)*.25
 # The public mass wrapper centers the sample before the ABI call.
 centered=dict(a,z=a['z']-kink)
 status,mass=core.run('mass_'+method,centered,{},threads=1);assert status==0
 options=dict(core_intervals=6,grid_points=101,lambda_points=5,refinement_points=3)
 fixed=dict(options,B=mass['B'],sigma=mass['sigma'],a=mass['a'] if method=='voigt' else 0.,kink=kink)
 status,raw=core.run('npmle',a,fixed,threads=1);assert status==0,(name,raw['status']);raw.pop('threads_used')
 case=dict(id='npmle_'+name,entry='npmle',options=options,mass_options=dict(method=method,kink=kink),columns=save_arrays(args.root,name,a),expected=encode(status,raw))
 cases.append(case);print(name,'PASS',flush=True)
write_json(args.root/'manifest.json.gz',dict(format='spectral-xlang-v1',cases=cases))
