"""Portable simulated inputs, raw binary format, and exact-output encoding."""
import hashlib, struct
from pathlib import Path
import numpy as np


def simulate(n, seed=1, shape="base", rd=False):
    # Counter-based SplitMix64, then a sum of 12 dyadic uniforms (no libm).
    # Every vector operation has its specified order; see README for format.
    counters = np.arange(n, dtype=np.uint64)
    def uniform(stream):
        x = counters + np.uint64((seed + stream * 0x9e3779b97f4a7c15) & ((1<<64)-1))
        x = (x ^ (x >> 30)) * np.uint64(0xbf58476d1ce4e5b9)
        x = (x ^ (x >> 27)) * np.uint64(0x94d049bb133111eb)
        x = x ^ (x >> 31)
        return (x >> 11).astype(np.float64) * (2.0**-53)
    def normal(stream):
        x = np.zeros(n)
        for j in range(12):
            x += uniform(stream + j)
        return x - 6.0
    latent = normal(1)
    noise = normal(20)
    if rd:
        d = (uniform(40) < np.where(latent >= 0, 0.8, 0.2)).astype(float)
        return {"w": latent + 0.12*noise, "d": d, "y": 1.5*d + 0.25*latent + normal(60)*0.2}
    member = uniform(40) < 0.2
    z = np.where(member, 0.0, latent) + 0.04*noise
    if shape == "skew":
        z = np.where(member, 0.0, latent*latent*0.4 - 0.35) + 0.08*noise
    elif shape == "refusal":
        z = latent
    return {"z": z, "x": np.where(member, 2.0, 0.5 + latent*0.2) + normal(60)*0.1}


def encode(status, raw):
    def number(x):
        return {"hex": float(x).hex(), "bits": struct.pack(">d", x).hex()}
    return {"status": status, "scalars": {k: number(v) for k,v in raw.items() if isinstance(v, (float,int))},
            "vectors": {k: [number(x) for x in v] for k,v in raw.items() if isinstance(v,np.ndarray)},
            "strings": {k:v for k,v in raw.items() if isinstance(v,str)}}


def save_arrays(root, name, arrays):
    columns = {}
    for key, value in arrays.items():
        path = root/(name+"."+key+".f64le")
        value.astype("<f8").tofile(path)
        columns[key] = {"file": path.name, "n": len(value), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return columns


def load_arrays(root, case):
    if "recipe" in case:
        a = simulate(**case["recipe"])
    else:
        a = {}
        for k, column in case["columns"].items():
            path = root/column["file"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == column["sha256"]
            a[k] = np.fromfile(path, dtype="<f8").astype(np.float64, copy=False)
            assert len(a[k]) == column["n"]
    return a


def public(case, arrays, threads=1):
    opts = dict(case["options"], threads=threads)
    entry = case["entry"]
    weights = arrays.get("weights")
    if entry == "rd":
        import specrd
        return specrd.rd(arrays["w"], arrays["d"], arrays["y"], weights=weights, **opts)
    import specbunch
    if entry == "npmle":
        fit=specbunch.mass(arrays["z"],weights=weights,threads=threads,**case["mass_options"])
        return specbunch.npmle(fit,arrays["z"],weights=weights,**opts)
    if entry == "bmean":
        return specbunch.buncher_mean(arrays["z"], arrays["x"], weights=weights, **opts)
    return specbunch.mass(arrays["z"], weights=weights, method=entry[5:], **opts)


def read_json(path):
    import gzip, json
    with gzip.open(path,"rt",encoding="utf-8") as f:
        return json.load(f)


def write_json(path, value):
    import gzip, json
    with gzip.open(path,"wt",encoding="utf-8",compresslevel=6) as f:
        json.dump(value,f,separators=(",",":"))
        f.write("\n")
