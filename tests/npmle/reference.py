"""Independent simulated-data reference: SciPy splines, likelihood and optimizer.

No package/core calls. Coefficients are nonnegative, with total smooth mass fixed.
Gaussian polynomial integrals are analytic. Voigt uses analytic Cauchy integrals
and high-order Gaussian expectation, independent of the core's adaptive method.
"""
import numpy as np
from scipy.interpolate import BSpline, PPoly
from scipy.special import ndtr, roots_hermitenorm
from scipy.optimize import minimize
from scipy.linalg import solve, qr, solve_triangular
from math import comb

class Side:
    def __init__(self,knots):
        self.knots=np.asarray(knots);self.t=np.r_[np.repeat(knots[0],3),knots,np.repeat(knots[-1],3)]
        self.m=len(self.t)-4
        self.area=(self.t[4:]-self.t[:-4])/4
        self.spl=BSpline(self.t,np.eye(self.m),3,extrapolate=False)
        self.anti=self.spl.antiderivative()
        self.poly=np.zeros((len(knots)-1,self.m,4))
        for j in range(self.m):
            pp=PPoly.from_spline((self.t,np.eye(self.m)[j],3))
            for i,lo in enumerate(knots[:-1]):self.poly[i,j]=pp.c[::-1,np.searchsorted(pp.x,lo,side='right')-1]/self.area[j]
    def latent(self,x):
        x=np.asarray(x);inside=(x>=0)&(x<=self.knots[-1]);p=np.zeros((len(x),self.m))
        p[inside]=self.spl(x[inside])/self.area
        c=self.anti(np.clip(x,0,self.knots[-1]))/self.area
        return p,c
    def gaussian(self,x,sigma):
        x=np.asarray(x)
        if sigma==0:return self.latent(x)
        pdf=np.zeros((len(x),self.m));cdf=pdf.copy()
        for i,(lo,hi) in enumerate(zip(self.knots[:-1],self.knots[1:])):
            h=hi-lo;coef=self.poly[i];active=np.flatnonzero(np.any(coef!=0,axis=1))
            done=x>hi+9*sigma;cdf[np.ix_(done,active)]+=sum(coef[active,k]*h**(k+1)/(k+1) for k in range(4))
            mask=(x>=lo-9*sigma)&(~done);y=x[mask];v0=(lo-y)/sigma;v1=(hi-y)/sigma
            p0=np.exp(-v0*v0/2)/np.sqrt(2*np.pi);p1=np.exp(-v1*v1/2)/np.sqrt(2*np.pi)
            moments=[np.where(v0>=0,ndtr(-v0)-ndtr(-v1),ndtr(v1)-ndtr(v0)),p0-p1]
            for k in range(2,5):moments.append((k-1)*moments[k-2]+v0**(k-1)*p0-v1**(k-1)*p1)
            e=[sum(comb(k,q)*(y-lo)**(k-q)*sigma**q*moments[q] for q in range(k+1)) for k in range(5)]
            pdf[np.ix_(mask,active)]+=sum(e[k][:,None]*coef[active,k] for k in range(4))
            cdf[np.ix_(mask,active)]+=sum((h**(k+1)*ndtr((y-hi)/sigma)+e[k+1])[:,None]*coef[active,k]/(k+1) for k in range(4))
        return np.maximum(pdf,0),np.clip(cdf,0,1)
    def cauchy(self,x,a):
        x=np.asarray(x);pdf=np.zeros((len(x),self.m));cdf=pdf.copy()
        for i,(lo,hi) in enumerate(zip(self.knots[:-1],self.knots[1:])):
            h=hi-lo;coef=self.poly[i];active=np.flatnonzero(np.any(coef!=0,axis=1));z=x-lo+1j*a
            far=np.abs(h/z)<.5;near=~far;mom=[]
            previous=np.zeros(len(x),complex)
            for k in range(5):
                value=np.zeros(len(x),complex)
                if k==0:value[near]=np.log(h-z[near])-np.log(-z[near])
                else:value[near]=h**k/k+z[near]*previous[near]
                ratio=h/z[far];term=np.ones(len(ratio),complex);series=term/(k+1)
                for q in range(1,50):term*=ratio;series+=term/(k+q+1)
                value[far]=-h**(k+1)/z[far]*series
                mom.append(value.imag/np.pi);previous=value
            pdf[:,active]+=sum(mom[k][:,None]*coef[active,k] for k in range(4))
            end=.5+np.arctan((x-hi)/a)/np.pi
            cdf[:,active]+=sum((h**(k+1)*end+mom[k+1])[:,None]*coef[active,k]/(k+1) for k in range(4))
        return np.maximum(pdf,0),np.clip(cdf,0,1)

class Basis:
    def __init__(self,left,right):
        self.left=Side(left);self.right=Side(right);self.area=np.r_[self.left.area,self.right.area];self.m=len(self.area)
    def latent(self,x):
        x=np.asarray(x);lp,lc=self.left.latent(-x);rp,rc=self.right.latent(x);lp[x>=0]=0;rp[x<0]=0
        return np.c_[lp,rp],np.c_[1-lc,rc]
    def gaussian(self,x,sigma):
        lp,lc=self.left.gaussian(-np.asarray(x),sigma);rp,rc=self.right.gaussian(x,sigma)
        return np.c_[lp,rp],np.c_[1-lc,rc]
    def cauchy(self,x,a):
        lp,lc=self.left.cauchy(-np.asarray(x),a);rp,rc=self.right.cauchy(x,a)
        return np.c_[lp,rp],np.c_[1-lc,rc]
    def observed(self,x,sigma,a,order=1024):
        x=np.asarray(x);pdf=np.zeros((len(x),self.m+1));cdf=pdf.copy()
        if a==0:
            pdf[:,:-1],cdf[:,:-1]=self.gaussian(x,sigma)
            pdf[:,-1]=np.exp(-.5*(x/sigma)**2)/(sigma*np.sqrt(2*np.pi));cdf[:,-1]=ndtr(x/sigma)
        elif sigma==0:
            pdf[:,:-1],cdf[:,:-1]=self.cauchy(x,a)
            pdf[:,-1]=a/(np.pi*(x*x+a*a));cdf[:,-1]=.5+np.arctan(x/a)/np.pi
        else:
            nodes,weights=roots_hermitenorm(order);weights/=np.sqrt(2*np.pi)
            keep=weights>1e-20;nodes=nodes[keep];weights=weights[keep]
            for start in range(0,len(x),32):
                stop=min(len(x),start+32);shift=(x[start:stop,None]-sigma*nodes).ravel()
                p,c=self.cauchy(shift,a);shape=(stop-start,len(nodes),self.m)
                pdf[start:stop,:-1]=np.einsum('ijk,j->ik',p.reshape(shape),weights)
                cdf[start:stop,:-1]=np.einsum('ijk,j->ik',c.reshape(shape),weights)
                pdf[start:stop,-1]=(a/(np.pi*(shift.reshape(stop-start,-1)**2+a*a)))@weights
                cdf[start:stop,-1]=(.5+np.arctan(shift.reshape(stop-start,-1)/a)/np.pi)@weights
        return pdf,cdf

def design(z,B,sigma,a,basis,width,weights=None,exact_threshold=1000,order=1024):
    z=np.asarray(z);w=np.ones(len(z)) if weights is None else np.asarray(weights);keep=np.isfinite(z)&np.isfinite(w)&(w>0);z=z[keep];w=w[keep]
    neff=w.sum()**2/(w@w)
    if len(z)<=exact_threshold:
        index=np.argsort(z,kind='stable');p,_=basis.observed(z[index],sigma,a,order);F=(1-B)*p[:,:-1];atom=B*p[:,-1];count=w[index]*neff/w.sum()
    else:
        key,inverse=np.unique(np.floor(z/width),return_inverse=True);count=np.bincount(inverse,weights=w)*neff/w.sum()
        if a>0:
            nodes,quadrature_weights=np.polynomial.legendre.leggauss(4)
            x=(key[:,None]+.5+.5*nodes)*width
            density,_=basis.observed(x.ravel(),sigma,a,order)
            delta=np.einsum('ijk,j->ik',density.reshape(len(key),4,basis.m+1),quadrature_weights)*width/2
        else:
            edges,which=np.unique(np.r_[key*width,(key+1)*width],return_inverse=True)
            _,c=basis.observed(edges,sigma,a,order);delta=np.maximum(0,c[which[len(key):]]-c[which[:len(key)]])
        F=(1-B)*delta[:,:-1];atom=B*delta[:,-1]
    P=np.zeros((basis.m,basis.m));off=0;factors=[]
    for side in (basis.left,basis.right):
        D=np.diff(np.eye(side.m),2,axis=0)*(1-B)/side.area
        P[off:off+side.m,off:off+side.m]=D.T@D;block=np.zeros((len(D),basis.m));block[:,off:off+side.m]=D;factors.append(block);off+=side.m
    return F,atom,count,P,neff,np.vstack(factors)

def fit_lambda(F,atom,count,P,neff,D,lam,initial):
    def objective(p):
        g=atom+F@p
        if np.any(g<=0):return 1e300,np.zeros(len(p))
        return (-count@np.log(g)+.5*lam*(p@P@p))/neff,(-F.T@(count/g)+lam*P@p)/neff
    result=minimize(objective,initial,jac=True,method='SLSQP',bounds=[(0,1)]*len(initial),constraints={'type':'eq','fun':lambda p:p.sum()-1,'jac':lambda p:np.ones(len(p))},options={'ftol':1e-13,'maxiter':2000})
    p=result.x
    # Equality-constrained Newton polish on the active face; add violated bounds.
    for _ in range(80):
        f,grad=objective(p);g=atom+F@p;H=((F.T*(count/g**2))@F+lam*P)/neff
        mean=p@grad;free=np.flatnonzero((p>1e-10)|(grad<mean-1e-9))
        A=np.block([[H[np.ix_(free,free)],np.ones((len(free),1))],[np.ones((1,len(free))),np.zeros((1,1))]])
        rhs=np.r_[-grad[free],0.]
        try:d=solve(A,rhs,assume_a='sym')[:-1]
        except np.linalg.LinAlgError:break
        step=np.zeros_like(p);step[free]=d;alpha=min(1.,min([-p[j]/step[j] for j in free if step[j]<0]+[1.]))
        gd=grad@step
        if gd>=-1e-15:break
        for _ in range(50):
            trial=p+alpha*step;fn,_=objective(trial)
            if fn<=f+1e-4*alpha*gd:break
            alpha*=.5
        if alpha<1e-14:break
        p=trial
    g=atom+F@p;ll=count@np.log(g);information=(F.T*(count/g**2))@F
    weighted=F*(np.sqrt(count)/g)[:,None]
    augmented=np.vstack([weighted,np.sqrt(lam)*D]);scale=np.linalg.norm(augmented,axis=0)
    _,factor=qr(augmented/scale,mode='economic')
    influence=solve_triangular(factor.T,(weighted/scale).T,lower=True)
    edf=np.sum(influence**2)
    _,grad=objective(p);mean=p@grad;kkt=max(np.max(mean-grad),np.sum(p*np.maximum(grad-mean,0)))
    return dict(p=p,lambda_=lam,loglik=ll,edf=edf,bic=-2*ll+edf*np.log(neff),kkt=kkt)

def select(model,lambda_min=.01,lambda_max=10000,lambda_points=15,refinement_points=9):
    F,atom,count,P,neff,D=model;p=np.full(F.shape[1],1/F.shape[1]);path=[]
    grid=np.unique(np.exp(np.linspace(np.log(lambda_min),np.log(lambda_max),lambda_points)))
    for lam in grid:
        f=fit_lambda(F,atom,count,P,neff,D,lam,p);path.append(f);p=f['p']
    coarse=list(path)
    step=max(np.log(10.)/4,(np.log(lambda_max)-np.log(lambda_min))/(lambda_points-1))
    lower,upper=min(1e-6,lambda_min),max(1e10,lambda_max)
    while True:
        winner=int(np.argmin([f['bic'] for f in coarse]))
        if winner==0 and coarse[0]['lambda_']>lower*(1+1e-12):
            endpoint=coarse[0];lam=max(lower,np.exp(np.log(endpoint['lambda_'])-step))
            f=fit_lambda(F,atom,count,P,neff,D,lam,endpoint['p']);coarse.insert(0,f);path.append(f)
        elif winner==len(coarse)-1 and coarse[-1]['lambda_']<upper*(1-1e-12):
            endpoint=coarse[-1];lam=min(upper,np.exp(np.log(endpoint['lambda_'])+step))
            f=fit_lambda(F,atom,count,P,neff,D,lam,endpoint['p']);coarse.append(f);path.append(f)
        else:
            break
    best=coarse[winner]
    if 0<winner<len(coarse)-1:
        p=coarse[winner-1]['p']
        for lam in np.exp(np.linspace(np.log(coarse[winner-1]['lambda_']),np.log(coarse[winner+1]['lambda_']),refinement_points)):
            f=fit_lambda(F,atom,count,P,neff,D,lam,p);p=f['p'];path.append(f)
            if f['bic']<best['bic']:best=f
    return best,path
