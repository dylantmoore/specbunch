import math
import numpy as np
import pytest
import specbunch as sb

def fit(sigma=.03,a=0,kink=0):
    raw=dict(B=.08,sigma=sigma,a=a,friction='voigt' if a else 'gaussian',window_lo=12.,window_hi=45.)
    r=sb.Result(0,raw,'mass');r.kink=kink;return r

def sample():
    rng=np.random.default_rng(99);n=12000;member=rng.uniform(size=n)<.08;latent=rng.normal(0,.5,n)
    z=np.where(member,0,latent)+.03*rng.normal(size=n)
    x=np.where(member,.8,.3+.2*latent)+.1*rng.normal(size=n)
    return z,x

@pytest.mark.parametrize('degree',[1,2])
def test_local_mean_numpy_reference(degree):
    z,x=sample();w=1+.25*(np.arange(len(z))%3);bw=.0025
    r=sb.nonbuncher_mean(z,x,fit(),mu_b=.8,se_mu_b=.01,weights=w,degree=degree)
    idx=np.floor(z/bw).astype(int)+40;use=(idx>=0)&(idx<80)
    W=np.bincount(idx[use],weights=w[use],minlength=80);W2=np.bincount(idx[use],weights=w[use]**2,minlength=80)
    mean=np.bincount(idx[use],weights=w[use]*x[use],minlength=80)/W
    var=np.bincount(idx[use],weights=w[use]*(x[use]-mean[idx[use]])**2,minlength=80)/(W-W2/W)*W2/W**2
    lo=np.arange(-40,40)*bw;hi=lo+bw;center=(lo+hi)/2
    cdf=lambda v:np.array([.5*math.erfc(-t/(.03*math.sqrt(2))) for t in v])
    atom=.08*(cdf(hi)-cdf(lo))/bw;Z=np.column_stack([np.ones(80),center,center**2,cdf(center)-.5])
    coef=np.linalg.lstsq(Z,W/(w.sum()*bw)-atom,rcond=None)[0]
    share=atom/(atom+np.maximum(Z@coef,1e-12));sel=abs(center)<.05
    X=(1-share[sel,None])*np.column_stack([center[sel]**p for p in range(degree+1)])
    weight=W[sel]*(1-abs(center[sel])/.05);y=mean[sel]-.8*share[sel]
    inv=np.linalg.inv(X.T@(weight[:,None]*X));beta=inv@(X.T@(weight*y))
    meat=X.T@((weight**2*var[sel])[:,None]*X);deriv=-inv@(X.T@(weight*share[sel]))
    se=math.sqrt((inv@meat@inv)[0,0]+deriv[0]**2*.01**2)
    assert not r.refused
    assert r.m_NB==pytest.approx(beta[0],abs=2e-11)
    assert r.se==pytest.approx(se,abs=2e-11)
    assert 0<r.pi_kink<1
    assert len(r.raw['local_fit.center'])==40

@pytest.mark.parametrize('sigma,a',[(.03,0),(.02,.01),(0,.03)])
def test_constant_characteristic_and_threads(sigma,a):
    z,_=sample();x=np.full(len(z),.4);f=fit(sigma,a)
    r=sb.nonbuncher_mean(z,x,f,mu_b=.4,threads=1)
    assert r.m_NB==pytest.approx(.4,abs=1e-12)
    other=sb.nonbuncher_mean(z,x,f,mu_b=.4,threads=4)
    for key,value in r.raw.items():
        if key=='threads_used':continue
        if isinstance(value,str):assert value==other.raw[key]
        else:np.testing.assert_array_equal(np.asarray(value).view(np.uint64),np.asarray(other.raw[key]).view(np.uint64))

def test_eti_direct_transform_and_rate_pairs():
    z,_=sample();f=fit();r=sb.eti.estimate(f,z,net_below=1.34,net_above=1.)
    t=np.arange(751)*.06;ph=np.exp(-.5*.03**2*t*t)
    C=np.cos(t[:,None]*z).mean(axis=1);S=np.sin(t[:,None]*z).mean(axis=1)
    weights=np.ones(201);weights[1:-1:2]=4;weights[2:-1:2]=2
    m=((weights*C[:201]/ph[:201]).sum()*.06/3-.08*12)/math.pi
    tt=t[200:];phi=ph[200:];wt=phi**2/tt**2;J=np.sum(wt*tt*S[200:]/phi)/wt.sum()
    assert r.shoulder_midpoint==pytest.approx(m,abs=1e-12)
    assert r.shoulder_jump==pytest.approx(J,abs=1e-12)
    assert r.saez.e==pytest.approx(.08/m/math.log(1.34),abs=1e-12)
    q=sb.eti.estimate(f,z,tax_below=-.34,tax_above=0)
    assert r.saez.e==q.saez.e
    assert r.continuous.d>=r.saez.d

def test_refusals_flags_and_bad_inputs():
    z,x=sample();f=fit()
    r=sb.nonbuncher_mean(np.ones(100),np.ones(100),f,mu_b=.8)
    assert r.refused and r.refusal_reason=='insufficient_local_variation'
    r=sb.nonbuncher_mean(z,x,f,mu_b=20.)
    assert not r.refused and 'outside_outcome_range' in r.flags
    for kw in [dict(degree=3),dict(bandwidth=0),dict(se_mu_b=-1),dict(mu_b=float('nan'))]:
        with pytest.raises((ValueError,sb.CoreError)):sb.nonbuncher_mean(z,x,f,**({'mu_b':.8}|kw))
    with pytest.raises(ValueError):sb.nonbuncher_mean(z,x,f,mu_b=.8,kink=1)
    with pytest.raises(sb.CoreError):sb.eti.estimate(f,z,net_below=1,net_above=1)
    with pytest.raises(ValueError):sb.eti.estimate(f,z,tax_below=0,tax_above=.2,net_below=1,net_above=.8)

def test_first_stage_fit_objects():
    z,x=sample()
    f=sb.mass(z,method='joint',plot=False,threads=1)
    m=sb.buncher_mean(z,x,threads=1)
    assert not f.refused and not m.refused
    r=sb.nonbuncher_mean(z,x,f,m)
    assert not r.refused and r.mu_b==m.mu and r.se_mu_b==m.se and r.B==f.B
