"""Spectral bunching estimators, computed entirely by the bundled C++ core."""
from ._transport import C, DP, HEADER, Library, Result, CoreError, arrays, pointer, integer

__version__ = "0.1.0"
POLICIES = {"tiered": 0, "first": 2}

class _Data(C.Structure):
    _fields_ = [("n", C.c_size_t), ("z", DP), ("w", DP)]
class _XY(C.Structure):
    _fields_ = [("n", C.c_size_t), ("z", DP), ("y", DP), ("w", DP)]
class _Mass(C.Structure):
    _fields_ = HEADER + [("threads", C.c_int32), ("check_all_windows", C.c_int32), ("seed", C.c_uint64),
                        ("distribution_checks", C.c_int32), ("policy", C.c_int32), ("lp_iterations", C.c_int64), ("lp_retry_iterations", C.c_int64)]
class _Mean(C.Structure):
    _fields_ = HEADER + [("threads", C.c_int32), ("reserved_band", C.c_int32), ("reserved_visibility", C.c_int32), ("R", C.c_int32), ("seed", C.c_uint64)]

_core = Library("specbunch", "sb")
for _name in ("mass_gaussian", "mass_voigt", "mass_gaussian_jump"):
    _core.bind(_name, [C.POINTER(_Data), C.POINTER(_Mass), C.POINTER(C.c_void_p)], C.c_int)
_core.bind("mass", [C.POINTER(_Data), C.POINTER(_Mass), C.c_int32, C.c_int32, C.POINTER(C.c_void_p)], C.c_int)
_core.bind("bmean", [C.POINTER(_XY), C.POINTER(_Mean), C.POINTER(C.c_void_p)], C.c_int)
_core.bind("opts_mass_init", [C.POINTER(_Mass)], None)
_core.bind("opts_bmean_init", [C.POINTER(_Mean)], None)

_core.bind("default_threads", [], C.c_int)

def default_threads():
    """Return the shared core default, including OMP_NUM_THREADS override."""
    return _core.fn("default_threads")()

def version():
    """Return package version, core version, and core ABI version."""
    return {"package": __version__, "core": _core.fn("version")().decode(), "abi": _core.fn("abi_version")()}

def mass(z, kink=0.0, weights=None, method="cosine", friction="gaussian", policy="tiered", threads=None,
         seed=1, check_all_windows=False, lp_iterations=5000, lp_retry_iterations=31250, npmle=False,
         plot=True):
    """Estimate mass by weighted nonlinear least squares in transform levels.

    method="cosine" fits the empirical cosine transform; method="joint" also
    fits sine and allows a density jump. friction is "gaussian" or "voigt".
    Fixed inverse-variance weights use the transform variance; the solver is
    Gauss-Newton/Newton. Return ranked windows and log-linear OLS/LAD companions.
    result.transform contains t, t_scaled, cosine, cosine_se, cosine_fit; joint
    fits add sine_scaled, sine_scaled_se and sine_fit. n_points describes the
    selected fitting window, separate from the 401-point display grid.

    plot=True (default) draws the empirical cosine transform with a +/-2 SE band,
    the fitted curve over the whole grid and the shaded selected window into
    result.figure (a second panel shows t*S(t) for joint fits); the x-axis is
    t*sigma_pilot, or raw t with plot="raw". plot=False or None skips it.
    Matplotlib is the optional extra specbunch[plot]; without it a warning is
    issued and figure is None. plt.show() is never called. result.plot(frequency)
    draws a new figure. With npmle=True the NPMLE figure is result.npmle.figure.

    policy is tiered (default) or first; checks flag coherent answers. The neutral
    bootstrap seed is 1. Distribution checks follow selection unless
    check_all_windows=True. LP budgets count simplex iterations, not time.
    Missing/nonfinite rows and nonpositive weights are removed together.
    The joint method does not support weights. See docs/REFERENCE.md for all fields.
    """
    if method not in ("cosine","joint"): raise ValueError("method must be cosine or joint")
    if friction not in ("gaussian","voigt"): raise ValueError("friction must be gaussian or voigt")
    if policy not in POLICIES: raise ValueError("policy must be tiered or first")
    if not isinstance(check_all_windows, bool): raise TypeError("check_all_windows must be boolean")
    input_z=z
    z,w=arrays(z,weights)
    if z is None: raise ValueError("z is required")
    z=z-float(kink) if kink!=0 else z
    o=_Mass();_core.fn("opts_mass_init")(C.byref(o))
    o.seed=integer(seed,64,True);o.threads=integer(0 if threads is None else threads)
    o.policy=POLICIES[policy];o.check_all_windows=2 if check_all_windows else 0
    o.lp_iterations=integer(lp_iterations,64);o.lp_retry_iterations=integer(lp_retry_iterations,64)
    if min(o.lp_iterations,o.lp_retry_iterations)<=0:raise ValueError("LP iteration budgets must be positive")
    if not isinstance(npmle,bool):raise TypeError("npmle must be boolean")
    from ._plot import transform_xscale
    xscale=transform_xscale(plot)
    d=_Data(len(z),pointer(z),pointer(w))
    result=Result(*_core.call("mass",d,o,int(method=="joint"),int(friction=="voigt")),"mass")
    result.kink=float(kink)
    if xscale is not None:_attach_figure(result,result.plot,xscale,"mass")
    if npmle:result.npmle=None if result.refused else globals()["npmle"](result,input_z,weights=w,threads=threads)
    return result

def buncher_mean(z,x,kink=0.0,weights=None,threads=None,seed=1,draws=499):
    """Estimate the cosine-ratio buncher mean and return all ranked bands.

    Uses the first coherent band passing its fit test, otherwise the first
    coherent band with poor_fit. Coherence requires a positive denominator
    and an estimate within the retained outcome range. draws must be >=19.
    """
    z,x,w=arrays(z,x,weights)
    if z is None or x is None:raise ValueError("z and x are required")
    z=z-float(kink) if kink!=0 else z
    o=_Mean();_core.fn("opts_bmean_init")(C.byref(o))
    o.seed=integer(seed,64,True);o.threads=integer(0 if threads is None else threads);o.R=integer(draws)
    d=_XY(len(z),pointer(z),pointer(x),pointer(w))
    result=Result(*_core.call("bmean",d,o),"bmean");result.kink=float(kink)
    return result

def _attach_figure(result,draw,arg,name):
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        import warnings
        warnings.warn(f"matplotlib is not installed, so {name} returned results without a plot; "
                      "install specbunch[plot] or pass plot=False", stacklevel=3)
    else:result.figure=draw(arg)

__all__=["default_threads","mass","buncher_mean","version","Result","CoreError","POLICIES"]

class _Npmle(C.Structure):
    _fields_=HEADER+[(k,C.c_int32) for k in ("threads","core_intervals","grid_points","exact_threshold","max_iterations","lambda_points","refinement_points","reserved")]+[(k,C.c_double) for k in ("B","sigma","a","kink","histogram_width","lambda_min","lambda_max","quadrature_tolerance","optimality_tolerance")]
_core.bind("opts_npmle_init",[C.POINTER(_Npmle)],None)
_core.bind("npmle",[C.POINTER(_Data),C.POINTER(_Npmle),DP,C.c_size_t,C.POINTER(C.c_void_p)],C.c_int)

def _result_table(raw,prefix):
    columns={k[len(prefix):]:v for k,v in raw.items() if k.startswith(prefix)}
    try:
        import pandas as pd
        return pd.DataFrame(columns)
    except ImportError:return columns

class NpmleResult:
    """Fixed-parameter spline fit with histogram shares, curves and an optional figure.

    ``figure`` holds the matplotlib Figure drawn by ``npmle`` (None when plot=False
    or matplotlib is missing); ``plot(kind)`` draws a new one.
    """
    def __init__(self,status_code,raw,xlabel=None):
        self.status_code=status_code;self.raw=raw;self.status=raw["status"]
        self.figure=None;self.xlabel=xlabel
        self.histogram=_result_table(raw,"histogram.")
        self.curves=_result_table(raw,"curves.")
        self.lambda_path=_result_table(raw,"lambda_path.")
        self.lambda_=raw["lambda"]
    def __getattr__(self,name):
        if name in self.raw:return self.raw[name]
        if "curves."+name in self.raw:return self.raw["curves."+name]
        raise AttributeError(name)
    def plot(self,kind="both",xlabel=None):
        """Return a new matplotlib Figure: kind "both" (default), "pdf" or "cdf".

        Density view: data histogram, fitted observed density (bunchers included) and
        fitted latent density excluding the atom. CDF view: empirical, fitted observed
        and fitted latent CDF (jump of B at the kink). Never calls plt.show().
        """
        from ._plot import draw
        return draw(self,kind,xlabel)
    def __repr__(self):return f"NpmleResult(status={self.status!r}, lambda_={self.lambda_:.6g}, edf={self.edf:.6g})"

def npmle(fit,z,weights=None,bins=None,*,threads=None,core_intervals=12,grid_points=1001,
          exact_threshold=1000,max_iterations=200,lambda_min=.01,lambda_max=10000,
          lambda_points=15,refinement_points=9,histogram_width=None,
          quadrature_tolerance=1e-9,optimality_tolerance=1e-7,plot="both"):
    """Fit a convex penalised spline density with mass and friction held fixed.

    Uses the complete usable sample, exact observations when n<=exact_threshold,
    otherwise a histogram of width max(sigma,a)/20. Separate clamped cubic
    splines allow a weakly identified density jump. BIC selects lambda, with an
    interior refinement. Returns histogram and curves tables.
    Main options
    ------------
    weights supplies positive observation weights, filtered jointly with z.
    bins supplies finite increasing output edges and does not restrict the fit.
    plot="both" (default), "pdf" or "cdf" draws the fitted latent and observed
    distributions against the data into result.figure (matplotlib, optional
    extra specbunch[plot]; a missing matplotlib gives a warning, not an error);
    plot=False or None skips it. The figure is never shown with plt.show().
    Advanced (solver settings)
    --------------------------
    threads, core_intervals, grid_points, exact_threshold, max_iterations,
    lambda_min, lambda_max, lambda_points, refinement_points, histogram_width,
    quadrature_tolerance and optimality_tolerance control numerical evaluation.
    See docs/REFERENCE.md for their defaults and integration diagnostics.
    """
    if not isinstance(fit,Result) or fit.kind!="mass" or fit.refused:
        raise ValueError("npmle requires a coherent mass result")
    from ._plot import plot_kind
    kind=plot_kind(plot);name=getattr(z,"name",None)
    z,w=arrays(z,weights)
    if z is None:raise ValueError("z is required")
    edges=arrays(bins)[0] if bins is not None else None
    o=_Npmle();_core.fn("opts_npmle_init")(C.byref(o))
    o.B=fit.B;o.sigma=fit.sigma;o.a=fit.raw.get("a",0) if fit.raw["friction"]=="voigt" else 0;o.kink=fit.kink
    values=locals()
    for k in ("core_intervals","grid_points","exact_threshold","max_iterations","lambda_points","refinement_points"):
        setattr(o,k,integer(values[k]))
    o.threads=integer(0 if threads is None else threads)
    for k in ("lambda_min","lambda_max","quadrature_tolerance","optimality_tolerance"):setattr(o,k,float(values[k]))
    o.histogram_width=0 if histogram_width is None else float(histogram_width)
    d=_Data(len(z),pointer(z),pointer(w))
    result=NpmleResult(*_core.call("npmle",d,o,pointer(edges),0 if edges is None else len(edges)),
                       xlabel=name if isinstance(name,str) and name else None)
    if kind is not None:_attach_figure(result,result.plot,kind,"npmle")
    return result

__all__ += ["npmle","NpmleResult"]

from ._postestimation import nonbuncher_mean, PostestimationResult
from . import eti
__all__ += ["nonbuncher_mean", "PostestimationResult", "eti"]
