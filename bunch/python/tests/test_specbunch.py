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

@pytest.mark.parametrize("method", ["gaussian", "voigt", "gaussian_jump"])
def test_results(sample, method):
    z, x = sample
    r = sb.mass(z, method=method, threads=1)
    assert isinstance(r.B, float) and isinstance(r.sigma, float)
    assert r.refused == (r.status != 0)
    assert isinstance(r.flags, list) and hasattr(r, "windows")
    assert "not a CI" in r.summary()
    assert r.se_scale == "log(B)"
    for key, col in {k:v for k,v in r.raw.items() if k.startswith("windows.")}.items():
        assert isinstance(col,np.ndarray) and col.flags.owndata
    if method == "gaussian":
        assert r.simple_ols and r.simple_lad
    old = r.raw["B"]
    sb.mass(z, method=method, threads=1)
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
    with pytest.raises(sb.CoreError): sb.mass(z,weights=np.ones(len(z)),method="gaussian_jump")
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
