"""Independent ABI oracle: no imports of specbunch/specrd or their transport."""
import ctypes as ct
import numpy as np

P = ct.POINTER(ct.c_double)
class MassData(ct.Structure):
    _fields_ = [("n", ct.c_size_t), ("z", P), ("w", P)]
class MeanData(ct.Structure):
    _fields_ = [("n", ct.c_size_t), ("z", P), ("y", P), ("w", P)]
class RDData(ct.Structure):
    _fields_ = [("n", ct.c_size_t), ("w", P), ("d", P), ("y", P), ("wt", P)]
class MassOpts(ct.Structure):
    _fields_ = [("struct_size", ct.c_uint32), ("abi_version", ct.c_uint32), ("threads", ct.c_int32),
                ("check_all_windows", ct.c_int32), ("seed", ct.c_uint64), ("distribution_checks", ct.c_int32),
                ("policy", ct.c_int32), ("lp_iterations", ct.c_int64), ("lp_retry_iterations", ct.c_int64)]
class MeanOpts(ct.Structure):
    _fields_ = [("struct_size", ct.c_uint32), ("abi_version", ct.c_uint32), ("threads", ct.c_int32),
                ("reserved_band", ct.c_int32), ("reserved_visibility", ct.c_int32), ("R", ct.c_int32), ("seed", ct.c_uint64)]
class RDOpts(ct.Structure):
    _fields_ = [("struct_size", ct.c_uint32), ("abi_version", ct.c_uint32), ("threads", ct.c_int32),
                ("reserved0", ct.c_int32), ("seed", ct.c_uint64)]
class NpmleOpts(ct.Structure):
    _fields_ = [("struct_size",ct.c_uint32),("abi_version",ct.c_uint32)]+[(k,ct.c_int32) for k in ("threads","core_intervals","grid_points","exact_threshold","max_iterations","lambda_points","refinement_points","reserved")]+[(k,ct.c_double) for k in ("B","sigma","a","kink","histogram_width","lambda_min","lambda_max","quadrature_tolerance","optimality_tolerance")]
POLICIES = {"tiered":0,"first":2}

class Direct:
    def __init__(self, path, rd=False):
        self.lib = ct.CDLL(str(path))
        self.prefix = "srd" if rd else "sb"
        self.rd = rd
        for name, args, restype in [
            ("result_free", [ct.c_void_p], None), ("result_status", [ct.c_void_p], ct.c_int),
            ("result_names", [ct.c_void_p, ct.c_int], ct.c_char_p),
            ("result_scalar", [ct.c_void_p, ct.c_char_p, P], ct.c_int),
            ("result_vector", [ct.c_void_p, ct.c_char_p, ct.POINTER(ct.c_size_t)], P),
            ("result_string", [ct.c_void_p, ct.c_char_p], ct.c_char_p),
            ("last_error", [], ct.c_char_p), ("build_info", [], ct.c_char_p),
        ]:
            f = self.f(name); f.argtypes, f.restype = args, restype
        entries = [("rd", RDData, RDOpts, "opts_init")] if rd else [
            ("mass_gaussian", MassData, MassOpts, "opts_mass_init"),
            ("mass_voigt", MassData, MassOpts, "opts_mass_init"),
            ("mass_gaussian_jump", MassData, MassOpts, "opts_mass_init"),
            ("bmean", MeanData, MeanOpts, "opts_bmean_init")]
        if not rd:
            self.f("opts_npmle_init").argtypes=[ct.POINTER(NpmleOpts)]
            self.f("npmle").argtypes=[ct.POINTER(MassData),ct.POINTER(NpmleOpts),P,ct.c_size_t,ct.POINTER(ct.c_void_p)]
            self.f("npmle").restype=ct.c_int
        for name, data, opts, init in entries:
            self.f(name).argtypes = [ct.POINTER(data), ct.POINTER(opts), ct.POINTER(ct.c_void_p)]
            self.f(name).restype = ct.c_int
            self.f(init).argtypes = [ct.POINTER(opts)]; self.f(init).restype = None

    def f(self, name):
        return getattr(self.lib, self.prefix + "_" + name)

    def run(self, entry, arrays, options, threads=1):
        # Fixture arrays are already contiguous float64. No package code used.
        a = arrays
        ptr = lambda k: a[k].ctypes.data_as(P) if k in a else None
        if entry == "rd":
            data = RDData(len(a["w"]), ptr("w"), ptr("d"), ptr("y"), ptr("weights"))
            opts = RDOpts(); self.f("opts_init")(ct.byref(opts))
        elif entry == "npmle":
            data = MassData(len(a["z"]),ptr("z"),ptr("weights"))
            opts=NpmleOpts();self.f("opts_npmle_init")(ct.byref(opts))
        elif entry == "bmean":
            data = MeanData(len(a["z"]), ptr("z"), ptr("x"), ptr("weights"))
            opts = MeanOpts(); self.f("opts_bmean_init")(ct.byref(opts))
        else:
            data = MassData(len(a["z"]), ptr("z"), ptr("weights"))
            opts = MassOpts(); self.f("opts_mass_init")(ct.byref(opts))
        opts.threads = threads
        for k, v in options.items():
            setattr(opts, k, POLICIES[v] if k == "policy" else v)
        h = ct.c_void_p()
        try:
            if entry == "npmle":
                edges=np.ascontiguousarray(options.get("bins",[]),dtype=np.float64)
                rc=self.f(entry)(ct.byref(data),ct.byref(opts),edges.ctypes.data_as(P),len(edges),ct.byref(h))
            else:
                rc = self.f(entry)(ct.byref(data), ct.byref(opts), ct.byref(h))
            if rc:
                raise RuntimeError(f"{rc}: {self.f('last_error')().decode()}")
            out = {}
            for kind in range(3):
                for name in self.f("result_names")(h, kind).decode().splitlines():
                    if kind == 0:
                        value = ct.c_double()
                        assert self.f("result_scalar")(h, name.encode(), ct.byref(value)) == 0
                        out[name] = value.value
                    elif kind == 1:
                        n = ct.c_size_t()
                        value = self.f("result_vector")(h, name.encode(), ct.byref(n))
                        out[name] = np.ctypeslib.as_array(value, shape=(n.value,)).copy() if n.value else np.empty(0)
                    else:
                        out[name] = self.f("result_string")(h, name.encode()).decode()
            return self.f("result_status")(h), out
        finally:
            if h.value:
                self.f("result_free")(h)
