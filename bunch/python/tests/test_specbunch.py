import numpy as np
import pytest
import specbunch as sb

@pytest.fixture
def sample():
    rng = np.random.default_rng(42)
    z = np.where(rng.random(1000)<0.3, 0, rng.normal(size=1000)) + rng.normal(0,0.04,1000)
    return z, 2 + rng.normal(size=1000)

def test_version():
    assert sb.version()["abi"] == 1
    assert sb.version()["package"] == "0.1.0"

@pytest.mark.parametrize("method,friction", [("cosine","gaussian"),("cosine","voigt"),("joint","gaussian"),("joint","voigt")])
def test_results(sample, method, friction):
    z, x = sample
    r = sb.mass(z, method=method, friction=friction, threads=1)
    assert isinstance(r.B, float) and isinstance(r.sigma, float)
    assert (r.method, r.friction) == (method, friction)
    table={k[10:]:v for k,v in r.raw.items() if k.startswith("transform.")}
    t=table["t"]
    assert len(t)==401 and np.all(np.diff(t)>0)
    c=np.cos(z[:,None]*t[None,:]);expected=c.mean(axis=0)
    np.testing.assert_allclose(table["cosine"],expected,atol=2e-13,rtol=2e-13)
    np.testing.assert_allclose(table["cosine_se"],np.sqrt(c.var(axis=0)/len(z)),atol=2e-13,rtol=2e-13)
    np.testing.assert_array_equal(table["t_scaled"],t*r.sigma_pilot)
    assert t[0] <= np.min(r.raw["windows.window_lo"])
    assert t[-1] >= np.max(r.raw["windows.window_hi"])
    if not r.refused:
        phi=np.exp(-r.raw.get("a",0)*t-r.sigma**2*t*t/2)
        np.testing.assert_allclose(table["cosine_fit"],r.B*phi,rtol=2e-14,atol=1e-15)
        selected=r.raw["windows.selected"]==1
        assert r.n_points==r.raw["windows.n_points"][selected][0]
    if method=="joint":
        v=np.sin(z[:,None]*t[None,:])
        np.testing.assert_allclose(table["sine_scaled"],t*v.mean(axis=0),atol=2e-12,rtol=2e-13)
        np.testing.assert_allclose(table["sine_scaled_se"],t*np.sqrt(v.var(axis=0)/len(z)),atol=2e-12,rtol=2e-13)
        if not r.refused:np.testing.assert_allclose(table["sine_fit"],r.jump*phi,rtol=2e-14,atol=1e-15)
    assert "Weighted nonlinear least squares" in r.estimator
    assert r.refused == (r.status != 0)
    assert isinstance(r.flags, list) and hasattr(r, "windows")
    assert "not a CI" in r.summary()
    assert r.se_scale == "log(B)"
    for key, col in {k:v for k,v in r.raw.items() if k.startswith("windows.")}.items():
        assert isinstance(col,np.ndarray) and col.flags.owndata
    if method == "cosine":
        assert r.loglinear_ols and r.loglinear_lad
    old = r.raw["B"]
    sb.mass(z, method=method, friction=friction, threads=1)
    assert r.raw["B"] == old or np.isnan(old)

@pytest.mark.parametrize("policy", list(sb.POLICIES))
def test_policies(sample, policy):
    r = sb.mass(sample[0], policy=policy, threads=1)
    assert r.policy == policy

def test_bmean(sample):
    z,x=sample
    r=sb.buncher_mean(z,x,weights=np.ones(len(z)),threads=1,draws=19)
    assert r.mu == r.mu or np.isnan(r.mu)
    assert "pointwise_median" in r.raw and "windows.mu" in r.raw
    assert r.raw["weighted"] == 1

@pytest.mark.parametrize("bad", [np.ones((2,1000)), np.ones((1000,1))])
def test_shape_rejected(bad):
    with pytest.raises(ValueError,match="one-dimensional"): sb.mass(bad)

def test_validation(sample):
    z,x=sample
    with pytest.raises(ValueError,match="length"): sb.mass(z,weights=np.ones(3))
    with pytest.raises(ValueError): sb.mass(z,method="unknown")
    with pytest.raises(ValueError): sb.mass(z,policy="unknown")
    with pytest.raises(TypeError): sb.mass(z,bogus=1)
    with pytest.raises(ValueError): sb.mass(z,seed=-1)
    with pytest.raises(TypeError): sb.mass(z,seed=1.5)
    assert sb.mass(z[:5]).reason == "insufficient_data"
    with pytest.raises(sb.CoreError): sb.mass(z,weights=np.ones(len(z)),method="joint")
    with pytest.raises(sb.CoreError): sb.buncher_mean(z,x,draws=1)

def test_input_conversion_and_missing(sample):
    z,x=sample
    r=sb.mass(z.tolist(),threads=1)
    s=sb.mass(np.r_[np.nan,z],threads=1)
    assert s.raw["n"] == len(z)
    assert r.B == s.B or (np.isnan(r.B) and np.isnan(s.B))
    shifted=sb.mass(z+2,kink=2,threads=1)
    assert np.isclose(r.B,shifted.B,rtol=1e-9,equal_nan=True)
    strided=np.column_stack([z,z]).ravel()[::2]
    assert sb.mass(strided,threads=1).B == r.B

def test_thread_defaults(sample, monkeypatch):
    serial = "openmp=off" in sb._core.fn("build_info")().decode()
    used = lambda n: 1 if serial else n
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    assert sb.buncher_mean(*sample, draws=19).threads_used == used(sb.default_threads())
    for t in (1, 3, 16):
        assert sb.buncher_mean(*sample, draws=19, threads=t).threads_used == used(t)
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    assert sb.default_threads() == 3
    for t in (None, 0, -1):
        assert sb.buncher_mean(*sample, draws=19, threads=t).threads_used == used(3)
    assert sb.buncher_mean(*sample, draws=19, threads=16).threads_used == used(16)


@pytest.mark.slow
def test_window_checks_do_not_change_selection(sample):
    a=sb.mass(sample[0],threads=1)
    b=sb.mass(sample[0],threads=1,check_all_windows=True)
    for k in ("B","sigma","se_log_B","se_log_sigma","rank","tier","flags"):
        assert np.array_equal(np.asarray(a.raw[k]).view(np.uint64),np.asarray(b.raw[k]).view(np.uint64))
    assert np.sum(a.raw["windows.selected"]) == (not a.refused)
    assert np.array_equal(a.raw["windows.rank"],np.arange(1,len(a.raw["windows.rank"])+1))
    for obsolete in ("strict","fallback","tiered_u"):
        with pytest.raises(ValueError):sb.mass(sample[0],policy=obsolete)


def test_api_axes_validation():
    for obsolete in ("gaussian", "voigt", "gaussian_jump"):
        with pytest.raises(ValueError, match="method"):
            sb.mass([], method=obsolete)
    with pytest.raises(ValueError, match="friction"):
        sb.mass([], friction="unknown")
    assert sb.mass([]).method == "cosine"
    assert sb.mass([]).friction == "gaussian"


def test_weighted_transform(sample):
    z=sample[0];w=1+(np.arange(len(z))%3)*.25
    r=sb.mass(z,weights=w,threads=1)
    t=r.raw["transform.t"];c=np.cos(z[:,None]*t[None,:]);mean=np.sum(w[:,None]*c,axis=0)/w.sum()
    se=np.sqrt(np.sum(w[:,None]**2*(c-mean)**2,axis=0))/w.sum()
    np.testing.assert_allclose(r.raw["transform.cosine"],mean,atol=2e-13,rtol=2e-13)
    np.testing.assert_allclose(r.raw["transform.cosine_se"],se,atol=2e-13,rtol=2e-13)
    assert sb.mass([]).raw["transform.t"].size==0


def test_joint_voigt_checks(sample):
    r=sb.mass(sample[0],method="joint",friction="voigt",check_all_windows=True,threads=1)
    assert r.method=="joint" and r.friction=="voigt"
    for flag in ("cdf_inconsistent","moments_inconsistent","no_valid_distribution","distribution_check_inconclusive","residual_cdf_inconsistent"):
        assert r.raw[flag]==0
        assert np.all(r.raw["windows."+flag]==-1)
    if not r.refused:
        assert 0<r.B<1 and r.a>=0 and r.s2>=0 and r.sigma_eq>0
        assert np.isclose(r.sigma_eq,r.a+np.sqrt(r.a*r.a+r.s2))
    with pytest.raises(sb.CoreError,match="weights"):
        sb.mass(sample[0],weights=np.ones(len(sample[0])),method="joint",friction="voigt")
