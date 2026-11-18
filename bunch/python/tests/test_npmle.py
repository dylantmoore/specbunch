import numpy as np
import pytest
import specbunch as sb

def sample(n=1500):
    rng=np.random.default_rng(37)
    return np.where(rng.random(n)<.2,0,rng.normal(size=n))+rng.normal(0,.04,n)

def test_curves_histogram_and_threads():
    z=sample();fit=sb.mass(z,threads=1)
    bins=np.linspace(z.min()-.1,z.max()+.1,51)
    result=sb.npmle(fit,z,bins=bins,threads=1,grid_points=101,lambda_points=5,refinement_points=3)
    assert result.status=='converged'
    assert np.all(result.coefficients>=0)
    assert abs(result.empirical_cdf[-1]-1)<1e-14
    for key in ['latent_cdf','observed_cdf','empirical_cdf']:assert np.min(np.diff(getattr(result,key)))>=-1e-12
    i=np.searchsorted(result.grid,0);assert abs(result.latent_cdf[i]-result.latent_cdf[i-1]-fit.B)<1e-12
    h=result.raw
    np.testing.assert_array_equal(h['histogram.fitted_share'],h['histogram.nonbuncher_share']+h['histogram.buncher_share'])
    assert abs(h['histogram.observed_share'].sum()-1)<1e-12
    for t in (2,4):
        other=sb.npmle(fit,z,bins=bins,threads=t,grid_points=101,lambda_points=5,refinement_points=3)
        for k,v in result.raw.items():
            if k=='threads_used':continue
            if isinstance(v,str):assert v==other.raw[k]
            else:assert np.array_equal(np.asarray(v).view(np.uint64),np.asarray(other.raw[k]).view(np.uint64))

def test_convenience_weights_and_full_sample():
    z=sample();fit=sb.mass(z,threads=1,npmle=True)
    assert isinstance(fit.npmle,sb.NpmleResult)
    weights=1+np.arange(len(z))%3
    edges=np.array([-.1,0,.1])
    result=sb.npmle(fit,z,weights=weights,bins=edges,threads=1,grid_points=51,lambda_points=3,refinement_points=3)
    assert result.n==len(z) and result.grid[0]==z.min() and result.grid[-1]==z.max()
    expected=np.histogram(z,edges,weights=weights)[0]/weights.sum()
    np.testing.assert_allclose(result.raw['histogram.observed_share'],expected,rtol=0,atol=2e-15)
    assert np.sum(expected)<1
    assert sb.mass(z[:5],npmle=True).npmle is None

def test_validation_and_exact_observations():
    z=sample(400);fit=sb.mass(z,threads=1)
    out=sb.npmle(fit,z,threads=1,grid_points=51,lambda_points=3,refinement_points=3)
    assert out.exact_observations==1 and out.likelihood_rows==400
    with pytest.raises(ValueError):sb.npmle(sb.mass(z[:5]),z)
    with pytest.raises(sb.CoreError):sb.npmle(fit,z,bins=[1,0])
    with pytest.raises(sb.CoreError):sb.npmle(fit,z,lambda_max=np.inf)


def test_nonzero_kink_curve_coordinates():
    z=sample()+.125;fit=sb.mass(z,kink=.125,threads=1)
    out=sb.npmle(fit,z,threads=1,grid_points=51,lambda_points=3,refinement_points=3)
    i=np.searchsorted(out.grid,.125)
    assert out.grid[i-1]==np.nextafter(.125,-np.inf)
    assert out.grid[i]==.125
    assert abs(out.latent_cdf[i]-out.latent_cdf[i-1]-fit.B)<1e-12
    np.testing.assert_allclose(out.empirical_cdf,np.searchsorted(np.sort(z),out.grid,side='right')/len(z),rtol=0,atol=1e-14)


def test_default_plot_and_opt_out(monkeypatch):
    pytest.importorskip("matplotlib")
    from matplotlib.figure import Figure
    z=sample(800);fit=sb.mass(z,threads=1)
    small=dict(threads=1,grid_points=101,lambda_points=3,refinement_points=3)
    out=sb.npmle(fit,z,**small)
    assert isinstance(out.figure,Figure) and len(out.figure.axes)==2
    pdf,cdf=out.figure.axes
    assert pdf.get_xlabel()=="z" and pdf.get_legend() is not None and cdf.get_legend() is not None
    assert any("hat B" in t.get_text() for t in pdf.texts)
    assert len(sb.npmle(fit,z,plot="cdf",**small).figure.axes)==1
    for off in (False,None,"none"):assert sb.npmle(fit,z,plot=off,**small).figure is None
    assert isinstance(out.plot("pdf"),Figure)
    with pytest.raises(ValueError):sb.npmle(fit,z,plot="hist",**small)
    monkeypatch.setitem(__import__("sys").modules,"matplotlib",None)
    with pytest.warns(UserWarning,match="matplotlib"):missing=sb.npmle(fit,z,**small)
    assert missing.figure is None and missing.status=="converged"
