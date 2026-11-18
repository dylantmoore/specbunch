"""C-ABI regression checks for upward/downward expansion and capped searches.

Run with --library pointing to an external build. Uses simulated data only.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'xlang'))
from raw_ctypes import Direct
from fixtures import simulate

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--library', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    core = Direct(a.library)
    z = simulate(600, seed=29)['z']
    receipts = []
    for name, scale, lo, hi, edge in [
        ('upward', 1., .01, 1., False),
        ('downward', 1., 1e7, 1e8, False),
        ('upper_cap', 1e6, .01, 1e4, True),
        ('lower_cap', 1e-6, .01, 1e4, True),
        ('adjacent_endpoints', 1., 1., np.nextafter(1., 2.), False),
    ]:
        status, r = core.run('npmle', {'z': z * scale}, dict(
            B=.2, sigma=.04 * scale, a=0., core_intervals=6, grid_points=25,
            lambda_points=5, refinement_points=3, lambda_min=lo, lambda_max=hi))
        assert status == 0, (name, r['status'])
        assert bool(r['lambda_at_grid_edge']) == edge, (name, r['lambda'])
        path = r['lambda_path.lambda']
        assert np.all(np.diff(path) > 0) and len(path) < 150
        if name == 'upward': assert r['lambda'] > hi
        if name == 'downward': assert r['lambda'] < lo
        if name == 'upper_cap': assert r['lambda'] == 1e10
        if name == 'lower_cap': assert r['lambda'] == 1e-6
        receipts.append(dict(case=name, lambda_=r['lambda'], edge=edge, points=len(path)))
    a.output.write_text(json.dumps(dict(result='PASS', cases=receipts), indent=2) + '\n')
    print('PASS', receipts)


if __name__ == "__main__":
    main()
