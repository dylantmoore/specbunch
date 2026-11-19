"""Bitwise local checks and a selected-answer cross-platform tolerance contract."""
import atexit
import gzip
import json
import math
import os
import platform
from pathlib import Path

EXACT_FIELDS = {'status', 'tier', 'flags', 'fit_converged', 'coherent_values',
                'cdf_inconsistent', 'moments_inconsistent', 'no_valid_distribution',
                'distribution_check_inconclusive', 'residual_cdf_inconsistent', 'poor_fit'}
# Endpoints identify the selected window and must agree exactly, as do decisions.
EXACT_FIELDS |= {'window_lo', 'window_hi', 'band_lo', 'band_hi', 't_lo', 't_hi'}
ESTIMATES = {'B', 'sigma', 'mu', 'tau', 'm_NB', 'saez.e', 'saez.d', 'continuous.e', 'continuous.d'}
STRINGS = {'reason', 'refusal_reason', 'flags', 'flags_str', 'status', 'saez.status', 'continuous.status'}


def close(x, y, rtol=1e-8, atol=1e-12):
    return x == y or (math.isnan(x) and math.isnan(y)) or (
        math.isfinite(x) and math.isfinite(y) and abs(x-y) <= max(atol, rtol*max(abs(x), abs(y))))


def differences(actual, expected):
    """Yield every differing field, including missing fields and vector elements."""
    if actual['status'] != expected['status']:
        yield 'status', actual['status'], expected['status'], True, None
    for group in ('strings', 'scalars', 'vectors'):
        aa, bb = actual[group], expected[group]
        for key in sorted(aa.keys() | bb.keys()):
            protected = (group == 'strings' and key in STRINGS) or (group == 'scalars' and key in EXACT_FIELDS | ESTIMATES)
            if key not in aa or key not in bb:
                yield group+'.'+key, aa.get(key), bb.get(key), protected, None
                continue
            if group == 'strings':
                if aa[key] != bb[key]:
                    yield group+'.'+key, aa[key], bb[key], protected, None
                continue
            xs, ys = ([aa[key]], [bb[key]]) if group == 'scalars' else (aa[key], bb[key])
            if len(xs) != len(ys):
                yield group+'.'+key+'.length', len(xs), len(ys), False, None
            for i in range(max(len(xs), len(ys))):
                x, y = xs[i] if i < len(xs) else None, ys[i] if i < len(ys) else None
                field = group+'.'+key+(f'[{i}]' if group == 'vectors' else '')
                if x is None or y is None:
                    yield field, x, y, protected, None
                elif x['bits'] != y['bits']:
                    u, v = float.fromhex(x['hex']), float.fromhex(y['hex'])
                    rel = abs(u-v)/max(abs(u), abs(v), 1e-300) if math.isfinite(u) and math.isfinite(v) else (0.0 if close(u,v) else math.inf)
                    flip = protected and not (close(u,v) if key in ESTIMATES else close(u,v,0,0))
                    yield field, x, y, flip, rel


def compare(actual, expected, *, exact=False, report=None, fixture=None):
    diffs = list(differences(actual, expected))
    count = 1+len(expected['strings'])+len(expected['scalars'])+sum(map(len, expected['vectors'].values()))
    if exact:
        assert actual == expected, ('bitwise mismatch', diffs[:3])
    elif report is not None:
        report.add(fixture, diffs)
    else:
        assert not any(d[3] for d in diffs), ('selected-answer flip', [d for d in diffs if d[3]][:3])
    return dict(compared=count, bitwise=max(0,count-len(diffs)), maxrel=max((d[4] for d in diffs if d[4] is not None), default=0.0))


class CrossReport:
    """Stream bounded-size compressed evidence; also finalize on an early exception."""
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.stream = gzip.open(self.directory/'differences.jsonl.gz', 'wt', encoding='utf-8')
        self.fixtures = set()
        self.flips = set()
        self.rows = {}
        self.finished = False
        self.complete = False
        atexit.register(self.finish)

    def add(self, fixture, diffs):
        # Thread invariance is checked separately, bitwise. Report each fixture once.
        if fixture in self.fixtures: return
        self.fixtures.add(fixture)
        for field, actual, expected, flip, rel in diffs:
            if flip: self.flips.add(fixture)
            curve = field.startswith('vectors.curves.')
            curve_exceeded = None
            if curve and isinstance(actual, dict) and isinstance(expected, dict):
                curve_exceeded = not close(float.fromhex(actual['hex']), float.fromhex(expected['hex']), 1e-6, 1e-8)
            def displayed(value):
                if isinstance(value, dict) and "hex" in value:
                    number = float.fromhex(value["hex"])
                    return dict(value, value=number if math.isfinite(number) else value["hex"])
                return value
            row = dict(fixture=fixture, field=field, actual=displayed(actual), expected=displayed(expected),
                       relative_difference=rel if rel is None or math.isfinite(rel) else 'infinity',
                       selected_answer_flip=flip, curve_tolerance_exceeded=curve_exceeded)
            self.stream.write(json.dumps(row, allow_nan=False)+'\n')
            key = (fixture, field.split('[')[0])
            n, maximum = self.rows.get(key, (0, 0.0))
            self.rows[key] = (n+1, max(maximum, rel or 0.0))
        self.stream.flush()

    def finish(self):
        if self.finished: return
        self.finished = True
        self.stream.close()
        summary = ['## Cross-platform reproducibility', '',
                   f'Actual platform: {platform.system()} {platform.machine()}; reference: Linux x86_64 CI fixtures.',
                   f'Fixtures compared: {len(self.fixtures)}. **Selected-answer flips: {len(self.flips)}**. Comparison complete: {self.complete}.',
                   'Same-platform thread comparisons remain bitwise. Selected estimates: relative 1e-8 / absolute 1e-12; decisions and selected endpoints exact.',
                   'Other differences are informational. NPMLE curves use relative 1e-6 / absolute 1e-8 for reporting only.',
                   'The artifact differences.jsonl.gz lists every difference (hex values preserve exact floating-point values).', '',
                   '| Fixture | Field | Differences | Max relative difference |', '|---|---|---:|---:|']
        for (fixture, field), (n, rel) in sorted(self.rows.items()):
            summary.append(f'| {fixture} | {field} | {n} | {rel:.6g} |')
        text = '\n'.join(summary)+'\n'
        (self.directory/'summary.md').write_text(text, encoding='utf-8')
        (self.directory/'counts.json').write_text(json.dumps(dict(complete=self.complete, fixtures=len(self.fixtures), selected_answer_flips=len(self.flips), flipped_fixtures=sorted(self.flips), differences=sum(n for n,r in self.rows.values())), indent=2)+'\n')
        if os.environ.get('GITHUB_STEP_SUMMARY'):
            # GitHub limits each step summary to 1 MiB; full detail stays in artifact.
            with open(os.environ['GITHUB_STEP_SUMMARY'], 'a', encoding='utf-8') as f:
                f.write(text if len(text.encode()) < 900000 else '\n'.join(summary[:8])+'\nFull field summary is in the artifact.\n')

    def assert_pass(self):
        self.complete = True
        self.finish()
        assert not self.flips, ('Selected-answer flips', sorted(self.flips))
