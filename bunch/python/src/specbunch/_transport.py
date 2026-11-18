"""Private ctypes transport. No estimation or selection is performed here."""
import ctypes as C
import operator
from pathlib import Path
import numpy as np

DP = C.POINTER(C.c_double)
HEADER = [("struct_size", C.c_uint32), ("abi_version", C.c_uint32)]

class CoreError(RuntimeError):
    """The core rejected a call; includes the C return code and last-error text."""
    def __init__(self, code, message):
        self.code = code
        super().__init__(f"core error {code}: {message}")

def integer(value, bits=32, unsigned=False):
    value = operator.index(value)
    lo, hi = (0, 2**bits) if unsigned else (-2**(bits-1), 2**(bits-1))
    if not lo <= value < hi:
        raise ValueError(f"integer outside {bits}-bit range")
    return value

def arrays(*values):
    out = []
    n = None
    for value in values:
        if value is None:
            out.append(None)
            continue
        a = np.asarray(value, dtype=np.float64)
        if a.ndim != 1:
            raise ValueError("inputs must be one-dimensional")
        if n is None:
            n = len(a)
        if len(a) != n:
            raise ValueError("input lengths must match")
        out.append(np.ascontiguousarray(a))
    return out

def pointer(a):
    return None if a is None else a.ctypes.data_as(DP)

class Library:
    def __init__(self, package, prefix):
        self.prefix = prefix
        root = Path(__file__).parent
        files = [p for p in root.iterdir() if p.name in
                 (f"lib{package}.dylib", f"lib{package}.so", f"{package}.dll", f"lib{package}.dll")]
        if len(files) != 1:
            raise ImportError(f"expected one bundled {package} library, found {files}")
        self.path = files[0]
        self.lib = C.CDLL(str(self.path))
        for name in ("version", "build_info", "last_error"):
            self.bind(name, [], C.c_char_p)
        self.bind("abi_version", [], C.c_uint32)
        self.bind("result_status", [C.c_void_p], C.c_int)
        self.bind("result_names", [C.c_void_p, C.c_int], C.c_char_p)
        self.bind("result_scalar", [C.c_void_p, C.c_char_p, DP], C.c_int)
        self.bind("result_vector", [C.c_void_p, C.c_char_p, C.POINTER(C.c_size_t)], DP)
        self.bind("result_string", [C.c_void_p, C.c_char_p], C.c_char_p)
        self.bind("result_free", [C.c_void_p], None)
        if self.fn("abi_version")() != 1:
            raise ImportError("unsupported core ABI (expected 1)")

    def fn(self, name):
        return getattr(self.lib, self.prefix + "_" + name)

    def bind(self, name, args, result):
        f = self.fn(name)
        f.argtypes, f.restype = args, result
        return f

    def call(self, name, data, opts, *extra):
        h = C.c_void_p()
        try:
            rc = self.fn(name)(C.byref(data), C.byref(opts), *extra, C.byref(h))
            if rc:
                raise CoreError(rc, self.fn("last_error")().decode())
            raw = {}
            for kind in range(3):
                names = self.fn("result_names")(h, kind).decode().splitlines()
                for key in names:
                    k = key.encode()
                    if kind == 0:
                        value = C.c_double()
                        rc = self.fn("result_scalar")(h, k, C.byref(value))
                        if rc:
                            raise CoreError(rc, f"cannot read scalar {key}")
                        raw[key] = value.value
                    elif kind == 1:
                        n = C.c_size_t()
                        p = self.fn("result_vector")(h, k, C.byref(n))
                        raw[key] = np.ctypeslib.as_array(p, shape=(n.value,)).copy() if n.value else np.empty(0)
                    else:
                        raw[key] = self.fn("result_string")(h, k).decode()
            return self.fn("result_status")(h), raw
        finally:
            if h.value:
                self.fn("result_free")(h)

class Result:
    """Owned copies of every named C output; no live C result handle."""
    def __init__(self, status, raw, kind):
        self.status, self.raw, self.kind = status, raw, kind
        self.refused = status != 0
        self.reason = raw.get("reason", "")
        self.flags = list(filter(None, raw.get("flags_str", "").split(",")))
        self.policy = raw.get("policy", "")
        self.tier = raw.get("tier")
        self.window = raw.get("win", (raw.get("window_lo", raw.get("band_lo", raw.get("t_lo"))), raw.get("window_hi", raw.get("band_hi", raw.get("t_hi")))))
        self.band = self.window
        self.se = raw.get("se", raw.get("se_log_B", float("nan")))
        self.se_scale = "log(B)" if "se_log_B" in raw else "estimate"
        columns={k[8:]:v for k,v in raw.items() if k.startswith("windows.")}
        checks={"cdf_inconsistent","moments_inconsistent","no_valid_distribution","distribution_check_inconclusive","residual_cdf_inconsistent","poor_fit"}
        for key in checks & columns.keys():
            if kind=="mass":columns[key]=np.array([{-1:"not_evaluated",0:"fail",1:"pass"}[int(x)] for x in columns[key]])
        if kind=="bmean" and "poor_fit" in columns:columns["flags"]=np.where(columns["poor_fit"]>0,"poor_fit","")
        try:
            import pandas as pd
            self.windows=pd.DataFrame(columns)
        except ImportError:self.windows=columns
        self.simple={prefix:{k[len(prefix)+1:]:v for k,v in raw.items() if k.startswith(prefix+".")} for prefix in ("loglinear_ols","loglinear_lad")}
        self.loglinear_ols=self.simple["loglinear_ols"];self.loglinear_lad=self.simple["loglinear_lad"]
        columns={k[10:]:v for k,v in raw.items() if k.startswith("transform.")}
        try:
            import pandas as pd
            self.transform=pd.DataFrame(columns)
        except ImportError:self.transform=columns
        self.diagnostics = raw.copy()

    def __getattr__(self, name):
        try:
            return self.raw[name]
        except KeyError:
            raise AttributeError(name) from None

    def summary(self):
        fields = ("B", "sigma", "se") if self.kind == "mass" else (("mu", "se", "pointwise_median") if self.kind == "bmean" else ("tau", "se"))
        rows = [f"{self.kind} | {'refused' if self.refused else 'estimate'} (status {self.status})", "field          value"]
        rows += [f"{k:<14} {getattr(self, k, float('nan')):.8g}" for k in fields]
        rows += [f"band/window    {self.window}", f"SE scale       {self.se_scale}; not a CI"]
        if self.reason:
            rows.append(f"reason         {self.reason}")
        if self.flags:
            rows.append("flags          " + ", ".join(self.flags))
        return "\n".join(rows)

    def __repr__(self):
        return self.summary()
