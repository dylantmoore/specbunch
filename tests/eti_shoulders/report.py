"""Aggregate the completed deterministic simulated-data run; no intervals."""
import argparse,csv,statistics,math,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('csv',type=Path);p.add_argument('log',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
rows=list(csv.DictReader(a.csv.open()));assert len(rows)==32*100*6
assert len({(r['cell'],r['rep']) for r in rows})==3200
for cell in range(32):
 group=[r for r in rows if int(r['cell'])==cell]
 assert {int(r['rep']) for r in group}==set(range(100))
 for rep in range(100):
  assert {(r['shoulders'],r['rule']) for r in group if int(r['rep'])==rep}=={(s,t) for s in ('transform','linear','cubic') for t in ('saez','continuous')}
wall=float(re.search(r'wall_seconds=([\d.]+)',a.log.read_text()).group(1))
def table(items,cell=False):
 lines=['| '+('n | σ | friction | jump | ' if cell else '')+'Rule | Shoulders | Finite / reps | Bias | RMSE | Median absolute error |','|'+('---:|---:|---|---|' if cell else '')+'---|---|---:|---:|---:|---:|']
 keys=sorted({((r['cell'],) if cell else ())+ (r['rule'],r['shoulders']) for r in items},key=lambda k:(int(k[0]),*k[1:]) if cell else k)
 for key in keys:
  group=[r for r in items if (((r['cell'],) if cell else ())+(r['rule'],r['shoulders']))==key]
  errors=[float(r['e'])-.3 for r in group if math.isfinite(float(r['e']))]
  bias=statistics.mean(errors) if errors else math.nan;rmse=math.sqrt(statistics.mean([x*x for x in errors])) if errors else math.nan;med=statistics.median([abs(x) for x in errors]) if errors else math.nan
  r=group[0];prefix=f"{r['n']} | {float(r['sigma']):.2f} | {'Voigt' if r['voigt']=='1' else 'Gaussian'} | {r['jump']} | " if cell else ''
  lines.append(f"| {prefix}{r['rule']} | {r['shoulders']} | {len(errors)}/{len(group)} | {bias:.5f} | {rmse:.5f} | {med:.5f} |")
 return '\n'.join(lines)
text=f'''# ETI shoulder Monte Carlo

Completed simulated-data experiment: 3,200 independent samples, 100 reps per cell;
32 cells = n ∈ {{200,000, 1,000,000}} × σ ∈ {{.02, .03, .05, .08}} ×
Gaussian/mild Voigt × no jump/jump. Same sample for all three shoulder methods.
True elasticity e=.3, log net-rate ratio q=.3. No SE, CI, or coverage calculation.
Final-run wall time: **{wall/60:.2f} minutes ({wall:.2f} seconds)** with eight
OpenMP workers across repetitions and one core thread per fit, on macOS arm64.

## Headline (equal replication weight across cells)

{table(rows)}

Errors are computed only for finite estimates; the denominator always counts all
replications. Missing C0 estimates from negative transform shoulders are shown
in the finite count, not silently removed from the denominator. This is a
**conditional shoulder experiment**: B and friction are set to their true values,
and the transform band is fixed at [2/σ, 3.5/σ]. It does not measure mass-estimator
error, window selection, or the full estimated-first-stage pipeline.

## Data-generating process and interpretation

Let d=e q=.09, L=.6, and R=.6 (no jump) or .9 (jump).
Set B=d(L+R)/2 and τ=(1−B)/(L+R). The latent post-policy distribution has an atom
B at zero and smooth densities L exp(z/τ) on z<0 and R exp(−z/τ) on z>0.
The no-policy density fills [0,d] by the line from L to R, retains the left tail,
and translates the right tail to begin at d. A homogeneous log-income shift of
d above the bunching interval produces this latent distribution. Thus the Saez
trapezoid is exact at the true shoulders. C0 uses a different completion of the
same gap; its error includes that approximation when L≠R, even at true shoulders.
With exact shoulders, the C0 elasticity in the jump design is
{3*(.09*(.6+.9)/2)/(.6+.9+math.sqrt(.6*.9))/.3:.8f}.

Observed z adds N(0,σ²) friction, plus independent Cauchy scale **a=.025σ** in
Voigt cells. This is mild Voigt friction, not a test of heavier Cauchy tails.
All sampled tails enter the likelihood; there is no truncation or winsorization.
Seeds are 20261119 + 100000*cell + rep, using C++ mt19937_64 and standard-library
normal draws. Reproducibility is deterministic on the same platform/build;
standard-library random distributions need not be identical across platforms.

## Numerical settings and limits

Both NPMLE bases use the default 12 core intervals per side and independent
clamped sides, fixed B/friction, and the existing second-difference height penalty.
To meet the runtime budget, likelihood bins are σ/2 (default σ/20), the initial
λ grid has 5 points (default 15), refinement 3 (default 9), the curve output grid
has 25 points, and Voigt quadrature tolerance is 1e−6 (default 1e−9).
The automatic λ extensions and caps are unchanged. These results apply to these
declared numerical settings, not an exhaustive comparison of tuning choices.
A four-interval pilot and a heavier-tail pilot were stopped during runtime
planning and are excluded from every headline and cell summary. The final
12-interval mild-Voigt pilot at n=1M, σ=.08 took 26.1 seconds for all methods;
this pilot is not an extra replication in the reported summaries.

The transform jump is a finite-band sine projection. Its noise can make one
shoulder negative even under equal true limits; Saez can remain finite while C0
is unavailable. NPMLE constrains both sides nonnegative but may introduce spline
and penalty bias. Neither method's shoulder diagnostic certifies identification.

## Fit diagnostics (one count per NPMLE fit)

'''
for basis in ('linear','cubic'):
 g=[r for r in rows if r['shoulders']==basis and r['rule']=='saez']
 text+=f"- {basis}: {sum(int(r['status'])!=0 for r in g)}/{len(g)} nonzero solver/quadrature statuses; {sum(float(r['edge'])!=0 for r in g)} grid-edge fits; {sum(float(r['weak'])!=0 for r in g)} weak-shoulder flags.\n"
text+='\n## Cell results\n\n'+table(rows,True)+'\n\nSource: `specbunch/tests/eti_shoulders/mc.cpp` and `report.py`. Only compact scalar records were persisted; generated samples were held in memory.\n'
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
print(table(rows));print('wall_seconds',wall)
