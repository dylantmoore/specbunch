#!/usr/bin/env python3
import os, tempfile
"""Exact public-e() parity, one Stata process, all fixture/thread combinations.
Uses the existing immutable spectral-xlang-v1 baseline. No baseline regeneration.
Generated large columns (~40 MB) and logs live only in the isolated build root.
"""
import argparse, hashlib, json, os, subprocess, time, sys, struct, math
from pathlib import Path
from fixtures import read_json, load_arrays
B=Path(os.environ.get("SPECPKG_BUILD", str(Path(tempfile.gettempdir())/"spectral-build"))) / 'stata'
p=argparse.ArgumentParser()
p.add_argument('--build-root',type=Path,default=B)
p.add_argument('--fixtures',type=Path,default=B.parent/'xlang/fixtures')
p.add_argument('--platform',choices=['macarm64','macintel64'],default='macarm64')
p.add_argument('--prepare-only',action='store_true')
p.add_argument('--compare-only',action='store_true',help='Check saved outputs from a completed real-Stata run')
p.add_argument('--limit',type=int)
a=p.parse_args(); B=a.build_root; out=B/'xlang'/a.platform; out.mkdir(parents=True,exist_ok=True)
manifest=a.fixtures/'manifest.json.gz'; cases=read_json(manifest)['cases']
if a.limit: cases=cases[:a.limit]
start=time.monotonic()
if not a.compare_only:
 plus=out/'plus'; personal=out/'personal'; plus.mkdir(exist_ok=True); personal.mkdir(exist_ok=True)
 lines=['clear all','set more off','set matsize 11000','set linesize 200',f'sysdir set PLUS "{plus}"',f'sysdir set PERSONAL "{personal}"', 'display c(stata_version) " " c(MP) " " c(machine_type) " " c(os)']
 for pkg in ['specbunch']:
  lines += [f'net install {pkg}, from("{B}/net/{pkg}") replace',f'program _{pkg}_plugin, plugin using("{pkg}_{a.platform}.plugin")']
 lines += r'''
 program define _xlang_export
 args filename
 mata: xlang_export(st_local("filename"))
 end
 mata:
 void xlang_export(string scalar path) {
 real scalar f,i,j,d
 string rowvector sn,bs,vn,tn
 string scalar k,m,s
 real colvector v
 f=fopen(path,"w")
 fput(f,"status"+char(9)+strofreal(st_numscalar("e(core_status)")))
 sn=tokens(st_global("e(scalar_names)")); bs=tokens(st_global("e(scalar_bits)"))
 for(i=1;i<=cols(sn);i++) {
 d=st_numscalar("e("+subinstr(sn[i],".","_")+")")
 s=missing(d) ? bs[i] : strtrim(strofreal(d,"%21x"))
 fput(f,"s"+char(9)+sn[i]+char(9)+s)
 }
 vn=tokens(st_global("e(vector_names)"))
 for(i=1;i<=cols(vn);i++) {
 k="vector"+strtrim(strofreal(i)); v=st_matrix("e("+k+")"); bs=tokens(st_global("e(bits_"+k+")"))
 s="v"+char(9)+vn[i]
 for(j=1;j<=rows(v);j++) s=s+char(9)+(missing(v[j]) ? bs[j] : strtrim(strofreal(v[j],"%21x")))
 fput(f,s)
 }
 tn=tokens(st_global("e(string_names)"))
 for(i=1;i<=cols(tn);i++) {
 k=st_global("e(str_"+subinstr(tn[i],".","_")+")"); s=""
 for(j=1;j<=strlen(k);j++) s=s+substr("0123456789abcdef",floor(ascii(substr(k,j,1))/16)+1,1)+substr("0123456789abcdef",mod(ascii(substr(k,j,1)),16)+1,1)
 fput(f,"t"+char(9)+tn[i]+char(9)+s)
 }
 fclose(f)
 }
 end
 '''.strip().splitlines()
 lines += r'''
 program define _npmle_export
 args filename
 mata: npmle_export(st_local("filename"))
 end
 mata:
 void npmle_export(string scalar path) {
 real scalar f,i,j,d
 string rowvector sn,bs,vn,tn
 string scalar k,m,s
 real colvector v
 f=fopen(path,"w")
 fput(f,"status"+char(9)+strofreal(st_numscalar("e(npmle_status_code)")))
 sn=tokens(st_global("e(npmle_scalar_names)")); bs=tokens(st_global("e(npmle_scalar_bits)"))
 for(i=1;i<=cols(sn);i++) {
 d=st_numscalar("e(npmle_"+subinstr(sn[i],".","_")+")")
 s=missing(d) ? bs[i] : strtrim(strofreal(d,"%21x"))
 fput(f,"s"+char(9)+sn[i]+char(9)+s)
 }
 vn=tokens(st_global("e(npmle_vector_names)"))
 for(i=1;i<=cols(vn);i++) {
 k="npmle_vector"+strtrim(strofreal(i)); v=st_matrix("e("+k+")"); bs=tokens(st_global("e(npmle_bits"+strtrim(strofreal(i))+")"))
 s="v"+char(9)+vn[i]
 for(j=1;j<=rows(v);j++) s=s+char(9)+(missing(v[j]) ? bs[j] : strtrim(strofreal(v[j],"%21x")))
 fput(f,s)
 }
 tn=tokens(st_global("e(npmle_string_names)"))
 for(i=1;i<=cols(tn);i++) {
 k=st_global("e(npmle_"+subinstr(tn[i],".","_")+")"); s=""
 for(j=1;j<=strlen(k);j++) s=s+substr("0123456789abcdef",floor(ascii(substr(k,j,1))/16)+1,1)+substr("0123456789abcdef",mod(ascii(substr(k,j,1)),16)+1,1)
 fput(f,"t"+char(9)+tn[i]+char(9)+s)
 }
 fclose(f)
 }
 end
 '''.strip().splitlines()
 paths={}; temporary=[]
 for c in cases:
  if 'recipe' in c:
   key=json.dumps(c['recipe'],sort_keys=True)
   if key not in paths:
    arrays=load_arrays(a.fixtures,c); entry={}
    for k,v in arrays.items():
     f=out/(('rd' if c['entry']=='rd' else 'bunch')+'_'+k+'.f64le'); v.astype('<f8').tofile(f); temporary.append(f); entry[k]=f
    paths[key]=entry
   cols=paths[key]; n=c['recipe']['n']
  else:
   cols={k:a.fixtures/v['file'] for k,v in c['columns'].items()}; n=next(iter(c['columns'].values()))['n']
   for k,f in cols.items(): assert hashlib.sha256(f.read_bytes()).hexdigest()==c['columns'][k]['sha256']
  keys=list(cols); pkg='specrd' if c['entry']=='rd' else 'specbunch'
  lines += ['clear',f'set obs {n}']+[f'generate double {k} = .' for k in keys]
  lines += [f'plugin call _{pkg}_plugin '+ ' '.join(keys)+', load '+' '.join(f'"{cols[k]}"' for k in keys)]
  for t in [1,2,4,8,0]:
   (out/f'{c["id"]}_{t}.tsv').unlink(missing_ok=True)
   weight=(' [fw=weights]' if pkg=='specrd' else ' [pw=weights]') if 'weights' in cols else ''
   opts=' '.join(f'{k}({v})' for k,v in c['options'].items())+(f' threads({t})' if t else '')
   if c['entry']=='rd': cmd='specrd y d w'
   elif c['entry']=='bmean': cmd='specbunch buncher z x'
   elif c['entry'] in ('nonbuncher_mean','eti'):
    f=c['fixed_mass']
    # Fixed first-stage oracle inputs isolate postestimation in each public wrapper.
    lines += ['capture program drop _fixed_mass','program define _fixed_mass, eclass','ereturn clear','ereturn local cmd "specbunch"','ereturn local method "cosine"','ereturn local friction "'+('voigt' if f['a']>0 else 'gaussian')+'"','ereturn scalar core_status=0']
    lines += [f'ereturn scalar {k}={v:.17g}' for k,v in f.items()]
    lines += ['end','_fixed_mass']
    aliases={'mu_b':'mub','se_mu_b':'semub','bin_width':'binwidth','net_below':'netbelow','net_above':'netabove','integration_intervals':'integrationintervals'}
    opts=' '.join(f'{aliases.get(k,k)}({v})' for k,v in c['options'].items())+(f' threads({t})' if t else '')
    cmd='specbunch '+c['entry']+' z'+(' x' if c['entry']=='nonbuncher_mean' else '')
   elif c['entry']=='npmle':
    mo=' '.join(f'{k}({v})' for k,v in c['mass_options'].items())
    lines += ['specbunch mass z'+weight+', '+mo+(f' threads({t})' if t else '')]
    aliases={'core_intervals':'coreintervals','grid_points':'gridpoints','exact_threshold':'exactthreshold','max_iterations':'maxiterations','lambda_min':'lambdamin','lambda_max':'lambdamax','lambda_points':'lambdapoints','refinement_points':'refinepoints','histogram_width':'histogramwidth','quadrature_tolerance':'quadraturetolerance','optimality_tolerance':'optimalitytolerance'}
    opts=' '.join(f'{aliases.get(k,k)}({v})' for k,v in c['options'].items())+(f' threads({t})' if t else '')
    cmd='specbunch npmle';weight=''
   else: cmd=f'specbunch mass z'; opts+=' method('+('joint' if c['entry'] in ('mass_gaussian_jump','mass_joint_voigt') else 'cosine')+') friction('+('voigt' if c['entry'] in ('mass_voigt','mass_joint_voigt') else 'gaussian')+')'
   lines += [cmd+weight+', '+opts, f'{"_npmle_export" if c["entry"]=="npmle" else "_xlang_export"} "{out}/{c["id"]}_{t}.tsv"',f'display "DONE {c["id"]} {t}"']
 lines += ['display "XLANG_STATA_PASS"','exit, clear']
 plan=out/'plan.do'; plan.write_text('\n'.join(lines)+'\n')
 (B/'xlang/plan.do').write_text(f'do "{plan}"\n')
 if a.prepare_only: print(plan); sys.exit()
 log=out/'plan.log'; log.unlink(missing_ok=True)
 command=[os.environ.get('STATA_EXECUTABLE','stata-mp'),'-b','do',str(plan)]
 if a.platform=='macintel64': command=['arch','-x86_64']+command
 start=time.monotonic()
 process=subprocess.Popen(command,cwd=out,env=dict(os.environ,TMPDIR=str(B/'tmp')))
 last=''
 while time.monotonic()-start<14400:
  text=log.read_text(errors='replace') if log.exists() else ''
  done=[x for x in text.splitlines() if x.startswith('DONE ')]
  if done and done[-1]!=last: last=done[-1]; print(last,flush=True)
  if '\nXLANG_STATA_PASS\n' in text: break
  if '\nend of do-file\n' in text and 'XLANG_STATA_PASS' not in text: raise RuntimeError(text[-2500:])
  time.sleep(1)
 else:
  process.terminate()
  raise TimeoutError('Stata did not finish; inspect '+str(log))
 process.wait(timeout=30)
else:
 assert '\nXLANG_STATA_PASS\n' in (out/'plan.log').read_text(errors='replace'), 'Saved Stata run did not complete'
 temporary=[]
def binary(s):
 if 'X' not in s:
  assert not math.isfinite(struct.unpack('>d',bytes.fromhex(s))[0]), ('finite output became Stata missing',s)
  return s
 mant,exponent=s.split('X')
 exponent_value=int(exponent,16)
 # Stata displays the stored exponent field for subnormals (zero minus 1023).
 if mant[1]=='0' and exponent_value==-1023:exponent_value=-1022
 d=float.fromhex(mant[0]+'0x'+mant[1:]+'p'+str(exponent_value))
 return struct.pack('>d',d).hex()
count=0
for c in cases:
 expected=c['expected']
 for t in [1,2,4,8,0]:
  got={'scalars':{},'vectors':{},'strings':{}}
  for line in (out/f'{c["id"]}_{t}.tsv').read_text().splitlines():
   fields=line.split('\t'); tag=fields[0]
   if tag=='status': got['status']=int(fields[1])
   elif tag=='s': got['scalars'][fields[1]]=binary(fields[2])
   elif tag=='v': got['vectors'][fields[1]]=[binary(x) for x in fields[2:]]
   elif tag=='t': got['strings'][fields[1]]=bytes.fromhex(fields[2] if len(fields)>2 else '').decode()
  want_threads=t or int(os.environ.get('SPECTRAL_EXPECT_DEFAULT','4'))
  assert got['scalars'].pop('threads_used') == struct.pack('>d',want_threads).hex()
  assert got['status']==expected['status'],(c['id'],t,'status'); count+=1
  for kind in ['scalars','vectors','strings']:
   assert got[kind].keys()==expected[kind].keys(),(c['id'],t,kind,'keys',got[kind].keys() ^ expected[kind].keys())
   for k,want in expected[kind].items():
    if kind=='scalars': want=want['bits']; count+=1
    elif kind=='vectors': want=[v['bits'] for v in want]; count+=len(want)
    else: count+=1
    assert got[kind][k]==want,(c['id'],t,kind,k,str(got[kind][k])[:180],str(want)[:180])
receipt={'fixtures':len(cases),'calls':len(cases)*5,'bitwise_comparisons':count,'threads':[1,2,4,8,0],'platform':a.platform,'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'status':'PASS','saved_output_comparison':a.compare_only}
(out/'result.json').write_text(json.dumps(receipt,indent=2)+'\n'); print(json.dumps(receipt),flush=True)
for f in temporary: f.unlink()
