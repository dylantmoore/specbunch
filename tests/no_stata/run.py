"""Drive the loaded plugin through an in-memory SPI host; compare all xlang fields."""
import argparse, ctypes as ct, sys, struct, json
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'xlang'))
from fixtures import read_json, load_arrays, encode
from raw_ctypes import POLICIES
from compare import compare, CrossReport
p=argparse.ArgumentParser()
p.add_argument('--host',type=Path,required=True);p.add_argument('--plugin',type=Path,required=True)
p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--family',choices=['bunch','rd'],required=True)
p.add_argument('--exact',action='store_true');p.add_argument('--output',type=Path)
p.add_argument('--case', help='Optional exact fixture id for a focused smoke test')
a=p.parse_args(); h=ct.CDLL(str(a.host.resolve()))
h.host_open.argtypes=[ct.c_char_p];h.host_open.restype=ct.c_int
h.host_data.argtypes=[ct.c_int,ct.c_int];h.host_column.argtypes=[ct.c_int,ct.POINTER(ct.c_double)]
h.host_call.argtypes=[ct.c_int,ct.POINTER(ct.c_char_p)];h.host_macro.argtypes=[ct.c_char_p];h.host_macro.restype=ct.c_char_p
h.host_matrix.argtypes=[ct.c_int,ct.c_int];h.host_value.argtypes=[ct.c_int];h.host_value.restype=ct.c_double
assert h.host_open(str(a.plugin.resolve()).encode())==0

def call(*args):
    argv=(ct.c_char_p*len(args))(*(str(x).encode() for x in args))
    assert h.host_call(len(args),argv)==0,args

def macro(k):return h.host_macro(('_'+k).encode()).decode()
def numbers(k):return [struct.unpack('>d',bytes.fromhex(b))[0] for b in macro(k).split()]
report = CrossReport((a.output.parent if a.output else a.fixtures.parent.parent)/'cross-os/stata') if not a.exact else None
total=dict(cases=0,calls=0,compared=0,bitwise=0,maxrel=0.0)
for case in read_json(a.fixtures/'manifest.json.gz')['cases']:
    if a.case and case['id']!=a.case:continue
    if (case['entry']=='rd') != (a.family=='rd'):continue
    arrays=load_arrays(a.fixtures,case)
    keys=['w','d','y'] if a.family=='rd' else (['z','x'] if case['entry'] in ('bmean','nonbuncher_mean') else ['z'])
    if 'weights' in arrays:keys+=['weights']
    h.host_data(len(arrays[keys[0]]),len(keys))
    for j,k in enumerate(keys):h.host_column(j,arrays[k].ctypes.data_as(ct.POINTER(ct.c_double)))
    first_actual=None
    for t in (1,2,4,8,0):
        o=case['options']; entry=case['entry'].removeprefix('mass_')
        if entry in ('nonbuncher_mean','eti'):
            f=case['fixed_mass']
            tail=[o.get('mu_b'),o.get('se_mu_b',0),o.get('bandwidth',.05),o.get('degree',1),o.get('bin_width',.0025)] if entry=='nonbuncher_mean' else [f['window_lo'],f['window_hi'],o['net_below'],o['net_above'],o.get('integration_intervals',200)]
            call(entry,f['B'],f['sigma'],f['a'],f['kink'],t,*tail)
        elif entry=='npmle':
            mo=case['mass_options'];kink=mo.get('kink',0)
            call('run',('joint_voigt' if mo['friction']=='voigt' else 'gaussian_jump') if mo['method']=='joint' else mo['friction'],kink,t,1,0,1,1,5000,31250,2,1,499)
            names=macro('scalar_names').split();h.host_matrix(1,len(names));call('fetch','scalars','out');mass=dict(zip(names,numbers('bits')));call('release')
            call('npmle',mass['B'],mass['sigma'],mass.get('a',0) if mo['friction']=='voigt' else 0,kink,t,o.get('core_intervals',12),o.get('grid_points',1001),o.get('exact_threshold',1000),o.get('max_iterations',200),o.get('lambda_min',.01),o.get('lambda_max',10000),o.get('lambda_points',15),o.get('refinement_points',9),o.get('histogram_width',0),o.get('quadrature_tolerance',1e-9),o.get('optimality_tolerance',1e-7),'_none')
        else:
            call('run',entry,0,t,o.get('seed',1),POLICIES[o.get('policy','tiered')], 1,1,5000,31250,2,1,499)
        raw={}; scalar_names=macro('scalar_names').split()
        h.host_matrix(1,len(scalar_names));call('fetch','scalars','out')
        for j,(k,v) in enumerate(zip(scalar_names,numbers('bits'))):
            stored=h.host_value(j)
            assert (stored>=8e307 if not np.isfinite(v) else stored==v),('matrix transport',k)
            raw[k]=v
        for j,k in enumerate(macro('vector_names').split(),1):
            n=int(macro('len'+str(j)));h.host_matrix(n,1);call('fetch',k,'out')
            raw[k]=np.asarray(numbers('bits'))
            assert len(raw[k])==n
            for i,v in enumerate(raw[k]):assert (h.host_value(i)>=8e307 if not np.isfinite(v) else h.host_value(i)==v)
        for j,k in enumerate(macro('string_names').split(),1):raw[k]=macro('str'+str(j))
        observed=raw.pop('threads_used')
        assert observed == t if t else observed >= 1
        actual=encode(int(macro('core_status')),raw)
        stats=compare(actual,case['expected'],exact=a.exact,report=report,fixture=case['id'])
        if first_actual is not None:compare(actual,first_actual,exact=True)
        else:first_actual=actual
        for k in ('compared','bitwise'):total[k]+=stats[k]
        total['maxrel']=max(total['maxrel'],stats['maxrel']);total['calls']+=1
        call('release')
    total['cases']+=1;print(case['id'],'PASS',flush=True)
assert total['cases'] > 0, 'No matching fixtures'
if report is not None: report.assert_pass()
total['result']='PASS';print(json.dumps(total))
if a.output:a.output.write_text(json.dumps(total,indent=2)+'\n')
