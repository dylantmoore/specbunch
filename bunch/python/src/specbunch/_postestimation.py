"""Transport for local nonbuncher means and shift-and-complete elasticity."""
from . import C, HEADER, _core, _Data, _XY, Result, arrays, pointer, integer, _result_table
from types import SimpleNamespace

class _Nonbuncher(C.Structure):
    _fields_=HEADER+[(k,C.c_int32) for k in ("threads","degree")]+[(k,C.c_double) for k in ("B","sigma","a","kink","mu_b","se_mu_b","bandwidth","bin_width")]
class _Eti(C.Structure):
    _fields_=HEADER+[(k,C.c_int32) for k in ("threads","integration_intervals")]+[(k,C.c_double) for k in ("B","sigma","a","kink","window_lo","window_hi","net_below","net_above")]
for name,opts,data in (("nonbuncher_mean",_Nonbuncher,_XY),("eti",_Eti,_Data)):
    _core.bind("opts_"+name+"_init",[C.POINTER(opts)],None)
    _core.bind(name,[C.POINTER(data),C.POINTER(opts),C.POINTER(C.c_void_p)],C.c_int)

class PostestimationResult:
    """Owned core outputs; SE is conditional and is not a confidence interval."""
    def __init__(self,status_code,raw):
        self.status_code=status_code;self.raw=raw;self.refused=status_code!=0
        self.flags=raw.get("flags","").split()
        self.local_fit=_result_table(raw,"local_fit.")
        for name in ("saez","continuous"):
            setattr(self,name,SimpleNamespace(**{k[len(name)+1:]:v for k,v in raw.items() if k.startswith(name+".")}))
    def __getattr__(self,name):
        if name in self.raw:return self.raw[name]
        raise AttributeError(name)
    def __repr__(self):return f"PostestimationResult(status={self.status!r}, flags={self.flags!r})"

def _mass_options(o,fit,kink):
    if not isinstance(fit,Result) or fit.kind!="mass" or fit.refused:
        raise ValueError("a coherent mass result is required")
    o.B=fit.B;o.sigma=fit.sigma;o.a=fit.raw.get("a",0) if fit.raw["friction"]=="voigt" else 0
    o.kink=fit.kink if kink is None else float(kink)
    if o.kink!=fit.kink:raise ValueError("kink must match the mass fit")

def nonbuncher_mean(z,x,mass_fit,buncher_mean_fit=None,*,mu_b=None,se_mu_b=0,
                    kink=None,weights=None,bandwidth=.05,degree=1,bin_width=.0025,threads=None):
    """Local nonbuncher mean at the kink, with the buncher mean held fixed.

    Pass a buncher_mean result or mu_b (and optionally se_mu_b). Triangular
    local-linear fitting is the default; degree=2 provides a quadratic check.
    The C++ core bins raw data and smooths the observed density. The result
    contains m_NB, se, pi_kink, flags and local_fit. SE is not a confidence
    interval: it conditions on mass/friction and density smoothing, propagating
    buncher-mean variance under independence. Report Gaussian and Voigt fits
    together because the mixture shares depend on the friction law.
    z is in log units, already transformed; kink defaults to mass_fit.kink.
    """
    if buncher_mean_fit is not None:
        if mu_b is not None:raise ValueError("supply buncher_mean_fit or mu_b, not both")
        if not isinstance(buncher_mean_fit,Result) or buncher_mean_fit.kind!="bmean" or buncher_mean_fit.refused:
            raise ValueError("a coherent buncher_mean result is required")
        mu_b=buncher_mean_fit.mu;se_mu_b=buncher_mean_fit.se
        if buncher_mean_fit.kink!=mass_fit.kink:raise ValueError("buncher mean and mass fit must have the same kink")
    if mu_b is None:raise ValueError("buncher_mean_fit or mu_b is required")
    z,x,w=arrays(z,x,weights)
    if z is None or x is None:raise ValueError("z and x are required")
    o=_Nonbuncher();_core.fn("opts_nonbuncher_mean_init")(C.byref(o));_mass_options(o,mass_fit,kink)
    o.mu_b=float(mu_b);o.se_mu_b=float(se_mu_b);o.bandwidth=float(bandwidth);o.degree=integer(degree);o.bin_width=float(bin_width);o.threads=integer(0 if threads is None else threads)
    return PostestimationResult(*_core.call("nonbuncher_mean",_XY(len(z),pointer(z),pointer(x),pointer(w)),o))

def _eti_estimate(fit,z,tax_below=None,tax_above=None,*,net_below=None,net_above=None,
                  kink=None,weights=None,integration_intervals=200,threads=None):
    """Saez trapezoid and continuous completion using CF-native shoulders.

    Supply both tax rates or both positive net-of-tax rates; the latter must
    decrease across the kink. z is log income, not income in levels. Returns
    saez.e/d, continuous.e/d, shoulder_left/right and flags. Shoulder estimates
    are finite-resolution CF functionals, not unrestricted pointwise densities.
    Negative shoulders flag the result and make continuous completion unavailable;
    Saez remains available when the shoulder midpoint is positive. No SE or CI
    is computed. See REFERENCE.md for structural assumptions.
    """
    if net_below is None and net_above is None:
        if tax_below is None or tax_above is None:raise ValueError("supply both tax rates or both net-of-tax rates")
        net_below=1-float(tax_below);net_above=1-float(tax_above)
    elif tax_below is not None or tax_above is not None or net_below is None or net_above is None:
        raise ValueError("supply exactly one complete rate pair")
    z,w=arrays(z,weights)
    if z is None:raise ValueError("z is required")
    o=_Eti();_core.fn("opts_eti_init")(C.byref(o));_mass_options(o,fit,kink)
    o.window_lo=fit.raw["window_lo"];o.window_hi=fit.raw["window_hi"]
    o.net_below=float(net_below);o.net_above=float(net_above);o.integration_intervals=integer(integration_intervals);o.threads=integer(0 if threads is None else threads)
    return PostestimationResult(*_core.call("eti",_Data(len(z),pointer(z),pointer(w)),o))
