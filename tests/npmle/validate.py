"""Full-size simulation-only validation against the independent SciPy reference."""
import argparse,json,time,sys,struct
from pathlib import Path
import numpy as np
import specbunch as sb
from reference import Basis,design,select
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--case',type=int);p.add_argument('--exclude-case',type=int);args=p.parse_args()
receipts=[]
cases=[(n,friction,jump) for n in (50000,1000000) for friction in ('gaussian','voigt') for jump in (False,True)]
for case,(n,friction,jump) in enumerate(cases):
 if args.case is not None and case!=args.case:continue
 if args.exclude_case is not None and case==args.exclude_case:continue
 start=time.monotonic();rng=np.random.default_rng(100+case)
 latent=np.abs(rng.normal(size=n))*np.where(rng.random(n)<(.65 if jump else .5),1,-1)
 latent[rng.random(n)<.2]=0
 z=latent+rng.normal(0,.04,n)
 if friction=='voigt':z+=.015*rng.standard_cauchy(n)
 method='gaussian_jump' if friction=='gaussian' and jump else friction
 mass=sb.mass(z,method=method,threads=4,seed=1)
 assert not mass.refused,(case,mass.reason)
 fit=sb.npmle(mass,z,threads=4)
 print(case,n,friction,jump,'native',fit.status,'lambda',fit.lambda_,'kkt',fit.kkt,flush=True)
 assert fit.status=='converged',(case,fit.status,fit.kkt)
 basis=Basis(fit.knots_left,fit.knots_right)
 model=design(z,fit.B,fit.sigma,fit.a,basis,fit.histogram_width)
 ref,path=select(model)
 errors={'lambda_relative':abs(ref['lambda_']/fit.lambda_-1),'loglik_absolute':abs(ref['loglik']-fit.loglik),'edf_absolute':abs(ref['edf']-fit.edf),'bic_absolute':abs(ref['bic']-fit.bic),'reference_kkt':ref['kkt']}
 print(case,'reference',errors,flush=True)
 assert errors['lambda_relative']<1e-10,errors
 assert errors['loglik_absolute']<.003 and errors['edf_absolute']<1e-4,errors
 edges=np.r_[fit.raw['histogram.bin_lo'],fit.raw['histogram.bin_hi'][-1]]
 _,cdf=basis.observed(edges,fit.sigma,fit.a)
 smooth=(1-fit.B)*(cdf[:,:-1]@ref['p']);atom=fit.B*cdf[:,-1]
 midpoint=(edges[:-1]+edges[1:])/2;latent_pdf,_=basis.latent(midpoint)
 histogram={'bin_lo':edges[:-1],'bin_hi':edges[1:],'observed_share':np.histogram(z,edges)[0]/n,'fitted_share':np.diff(smooth+atom),'nonbuncher_share':np.diff(smooth),'buncher_share':np.diff(atom),'latent_counterfactual_density':(1-fit.B)*(latent_pdf@ref['p'])}
 errors['table_max_absolute']=max(float(np.max(np.abs(v-fit.raw['histogram.'+k]))) for k,v in histogram.items())
 grid=fit.grid;lp,lc=basis.latent(grid);op,oc=basis.observed(grid,fit.sigma,fit.a)
 curves={'grid':grid,'latent_pdf':(1-fit.B)*(lp@ref['p']),'latent_cdf':(1-fit.B)*(lc@ref['p'])+fit.B*(grid>=0),'observed_pdf':(1-fit.B)*(op[:,:-1]@ref['p'])+fit.B*op[:,-1],'observed_cdf':(1-fit.B)*(oc[:,:-1]@ref['p'])+fit.B*oc[:,-1],'empirical_cdf':np.searchsorted(np.sort(z),grid,side='right')/n}
 errors['curves_max_absolute']=max(float(np.max(np.abs(v-fit.raw['curves.'+k]))) for k,v in curves.items())
 assert errors['table_max_absolute']<1e-6 and errors['curves_max_absolute']<2e-6,errors
 for threads in (1,2,8):
  other=sb.npmle(mass,z,threads=threads)
  assert other.raw.keys()==fit.raw.keys()
  for k,v in fit.raw.items():
   if k=='threads_used':continue
   u=other.raw[k]
   if isinstance(v,str):assert u==v,(case,threads,k)
   elif isinstance(v,np.ndarray):assert np.array_equal(v.view(np.uint64),u.view(np.uint64)),(case,threads,k)
   else:assert struct.pack('d',v)==struct.pack('d',u),(case,threads,k)
 receipt=dict(case=case,n=n,friction=friction,jump=jump,method=method,B=fit.B,sigma=fit.sigma,a=fit.a,lambda_=fit.lambda_,lambda_at_grid_edge=fit.lambda_at_grid_edge,errors=errors,thread_bitwise='PASS',status='PASS',seconds=time.monotonic()-start)
 receipts.append(receipt);args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(receipts,indent=2)+'\n');print('PASS',receipt,flush=True)
