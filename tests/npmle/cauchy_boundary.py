import argparse,json,numpy as np,specbunch as sb
from pathlib import Path
from reference import Basis,design,select
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
rng=np.random.default_rng(300);z=np.where(rng.random(5000)<.2,0,rng.normal(size=5000))+.04*rng.standard_cauchy(5000)
f=sb.Result(0,dict(B=.2,sigma=0.,a=.04,method='voigt'),'mass');f.kink=0
r=sb.npmle(f,z,threads=4,core_intervals=6,grid_points=101,lambda_points=5,refinement_points=3)
b=Basis(r.knots_left,r.knots_right);ref,_=select(design(z,r.B,r.sigma,r.a,b,r.histogram_width),lambda_points=5,refinement_points=3)
p,c=b.observed(r.grid,0,.04);err=max(np.max(abs((1-r.B)*(p[:,:-1]@ref['p'])+r.B*p[:,-1]-r.observed_pdf)),np.max(abs((1-r.B)*(c[:,:-1]@ref['p'])+r.B*c[:,-1]-r.observed_cdf)))
assert r.status=='converged';assert abs(ref['lambda_']/r.lambda_-1)<1e-10;assert abs(ref['loglik']-r.loglik)<.003;assert abs(ref['edf']-r.edf)<1e-4;assert err<2e-6
result=dict(status='PASS',lambda_=r.lambda_,loglik_error=abs(ref['loglik']-r.loglik),edf_error=abs(ref['edf']-r.edf),curve_error=err)
print(result);args.output.write_text(json.dumps(result,indent=2)+'\n')
