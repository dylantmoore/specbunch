"""Compact simulated fixtures with fixed first-stage inputs, from the direct C ABI."""
import argparse,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'xlang'))
from raw_ctypes import Direct
from fixtures import save_arrays,encode,write_json
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
a.root.mkdir(parents=True,exist_ok=True);lib=Direct(a.library);cases=[]
for law,sigma,scale in [('gaussian',.03,0.),('voigt',.02,.01),('cauchy',0.,.03)]:
 rng=np.random.default_rng(314);n=12000;member=rng.uniform(size=n)<.08;latent=rng.normal(0,.5,n)
 z=np.where(member,0,latent)+sigma*rng.normal(size=n)+scale*rng.standard_cauchy(n)
 x=np.where(member,.8,.3+.2*latent)+.1*rng.normal(size=n)
 arrays={'z':np.ascontiguousarray(z),'x':np.ascontiguousarray(x)}
 mass=dict(B=.08,sigma=sigma,a=scale,kink=0.,window_lo=12.,window_hi=45.)
 for degree in (1,2):
  if degree==2:arrays['weights']=np.ascontiguousarray(1+.25*(np.arange(n)%3))
  options=dict(mu_b=.8,se_mu_b=.01,degree=degree)
  ident=f'nonbuncher_{law}_{degree}';cols=save_arrays(a.root,ident,arrays)
  native={k:v for k,v in mass.items() if not k.startswith('window_')}|options
  status,raw=lib.run('nonbuncher_mean',arrays,native);raw.pop('threads_used')
  cases.append(dict(id=ident,entry='nonbuncher_mean',fixed_mass=mass,options=options,columns=cols,expected=encode(status,raw)))
 # Keep a weighted ETI fixture for each friction law.
 options=dict(net_below=1.34,net_above=1.)
 status,raw=lib.run('eti',arrays,mass|options);raw.pop('threads_used')
 cases.append(dict(id='eti_'+law,entry='eti',fixed_mass=mass,options=options,columns=cols,expected=encode(status,raw)))
# Refusals and a flagged answer are explicit fixtures.
for entry,label,z,options,B in [
 ('nonbuncher_mean','empty_local',np.full(100,2.),dict(mu_b=.8),.08),
 ('nonbuncher_mean','outside_range',np.linspace(-.2,.2,1000),dict(mu_b=20.),.08),
 ('eti','nonpositive_density',np.linspace(-2,2,1000),dict(net_below=1.34,net_above=1.),.8)]:
 arrays={'z':z,'x':np.linspace(0,1,len(z))};mass=dict(B=B,sigma=.03,a=0.,kink=0.,window_lo=12.,window_hi=45.)
 native=(mass if entry=='eti' else {k:v for k,v in mass.items() if not k.startswith('window_')})|options
 status,raw=lib.run(entry,arrays,native);raw.pop('threads_used');ident=entry+'_'+label
 cases.append(dict(id=ident,entry=entry,fixed_mass=mass,options=options,columns=save_arrays(a.root,ident,arrays),expected=encode(status,raw)))
write_json(a.root/'postestimation.json.gz',dict(format='spectral-xlang-v1',cases=cases))
print('Saved',len(cases),'postestimation fixtures')
