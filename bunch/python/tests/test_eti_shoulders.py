import numpy as np
import pytest
import specbunch as sb

def fixed():
    return sb.Result(0,dict(B=.08,sigma=.03,a=0.,kink=0.,friction='gaussian',window_lo=12.,window_hi=45.),'mass')

def test_npmle_eti_automatic_and_diagnostics():
    z=(np.arange(400,dtype=float)-200)/256
    import csv
    from pathlib import Path
    expected={row["basis"]:row for row in csv.DictReader((Path(__file__).parent/"eti_shoulders.csv").open())}
    f=fixed()
    for basis in ('linear','cubic'):
        n=sb.npmle(f,z,basis=basis,plot=False,core_intervals=4,grid_points=25,
                   lambda_points=3,refinement_points=3,threads=1,net_below=1.34,net_above=1)
        e=n.eti
        for key,actual in (("left",n.latent_left),("right",n.latent_right),("saez",e.npmle.saez.e),("continuous",e.npmle.continuous.e)):
            assert abs(actual-float(expected[basis][key]))<=max(1e-12,1e-8*abs(actual))
        assert e.shoulder_source=='transform' and e.npmle.shoulder_source=='npmle'
        assert e.npmle.basis==basis
        assert e.npmle.shoulder_left==n.latent_left==n.coefficients[0]
        assert e.npmle.saez.e==2*f.B/(n.latent_left+n.latent_right)/np.log(1.34)
        assert e.npmle.continuous.e==3*f.B/(n.latent_left+n.latent_right+np.sqrt(n.latent_left*n.latent_right))/np.log(1.34)
        path=n.lambda_path
        ll=np.asarray(path['latent_left']);rr=np.asarray(path['latent_right'])
        eligible=(np.asarray(path['converged'])!=0)&(np.asarray(path['bic'])<=n.bic+2)
        assert eligible.sum()==n.shoulder_lambda_count
        floor=max(1e-12,1e-6*(n.latent_left+n.latent_right)/2)
        assert n.shoulder_sensitivity_left==max(abs(ll[eligible]-n.latent_left)/max(n.latent_left,floor),default=0)
        f.npmle=n
        attached=sb.eti.estimate(f,z,net_below=1.34,net_above=1,threads=1)
        assert attached.npmle.saez.e==e.npmle.saez.e
        other=fixed();other.raw['B']=.09
        with pytest.raises(ValueError,match='match'):
            sb.eti.estimate(other,z,net_below=1.34,net_above=1,npmle=n,threads=1)
    with pytest.raises(ValueError):sb.npmle(f,z,basis='quadratic',plot=False)
