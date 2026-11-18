"""Independent NumPy/SciPy check of joint Voigt selected-window WNLS, simulated data only."""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0,str(Path(__file__).resolve().parent/'xlang'))
from raw_ctypes import Direct
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
core=Direct(args.library);rows=[]
for label,gaussian,cauchy,jump in [('gaussian',.04,0.,False),('voigt',.03,.02,False),('voigt_jump',.03,.02,True),('cauchy_jump',0.,.03,True)]:
    rng=np.random.default_rng(178);n=4000
    member=rng.random(n)<.25
    latent=rng.normal(size=n) if not jump else np.where(rng.random(n)<.7,1.,-1.)*rng.exponential(size=n)
    z=np.where(member,0.,latent)+gaussian*rng.normal(size=n)+cauchy*rng.standard_cauchy(n)
    start=time.monotonic();status,r=core.run('mass_joint_voigt',{'z':np.ascontiguousarray(z)},{'seed':1},threads=4);elapsed=time.monotonic()-start
    assert status==0,(label,r.get('reason'))
    t=np.linspace(r['window_lo'],r['window_hi'],int(r['n_points']));u=t*r['sigma_pilot']
    C=np.cos(z[:,None]*t);S=np.sin(z[:,None]*t);cos=C.mean(axis=0);sin=S.mean(axis=0)
    vc=(C*C).mean(axis=0)-cos*cos;vs=(S*S).mean(axis=0)-sin*sin;cs=(C*S).mean(axis=0)-cos*sin
    kap=cs/vs;conditional=vc-cs*cs/vs
    w1=np.sqrt(np.where(conditional>1e-12,n/conditional,n));w2=np.sqrt(np.where(vs>1e-12,n/vs,n))
    def model(x):
        B,J,a,b=x;phi=np.exp(-a*u-b*u*u/2);mc=B*phi;ms=J*phi/u
        return phi,mc,ms
    def residual(x):
        _,mc,ms=model(x);return np.r_[w1*(cos-kap*sin-mc+kap*ms),w2*(sin-ms)]
    def jacobian(x):
        phi,mc,ms=model(x);m1=mc-kap*ms
        d1=np.c_[phi,-kap*phi/u,-u*m1,-u*u*m1/2]
        d2=np.c_[np.zeros(len(u)),phi/u,-u*ms,-u*u*ms/2]
        return -np.r_[w1[:,None]*d1,w2[:,None]*d2]
    x=np.array([r['B'],r['jump']*r['sigma_pilot'],r['a']/r['sigma_pilot'],r['s2']/r['sigma_pilot']**2])
    fits=[least_squares(residual,start,jac=jacobian,bounds=([-np.inf,-np.inf,0,0],[np.inf]*4),ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=1000) for start in (x,np.array([.2,0,.3,.7]),np.array([.2,0,0.,1.]),np.array([.2,0,1.,0.]))]
    ref=min(fits,key=lambda f:np.dot(f.fun,f.fun));loss=np.dot(residual(x),residual(x));ref_loss=np.dot(ref.fun,ref.fun)
    assert loss<=ref_loss+1e-6*max(1.,ref_loss),(label,'objective',loss,ref_loss,x,ref.x)
    active=[0,1]+([2] if x[2]>0 else [])+([3] if x[3]>0 else [])
    jac=jacobian(x)[:,active];cov=np.linalg.inv(jac.T@jac);se=np.sqrt(cov[0,0])/x[0]
    assert np.isclose(r['se_log_B'],se,rtol=1e-4,atol=1e-8),(label,'standard error',r['se_log_B'],se)
    row=dict(case=label,B=r['B'],a=r['a'],s2=r['s2'],jump=r['jump'],core_objective=float(loss),reference_objective=float(ref_loss),max_parameter_difference=float(np.max(np.abs(x-ref.x))),se_log_B=r['se_log_B'],reference_se_log_B=float(se),seconds=elapsed,status='PASS')
    rows.append(row);print(json.dumps(row),flush=True)
args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps({'status':'PASS','objective_tolerance':'1e-6 * max(1, reference objective)','standard_error_tolerance':'relative 1e-4, absolute 1e-8','cases':rows},indent=2)+'\n')
