#!/usr/bin/env python3
import os, tempfile
"""Audit both Stata plugin families: symbols, private runtime, all Mach-O inputs."""
import hashlib,json,subprocess
from pathlib import Path
B=Path(os.environ.get("SPECPKG_BUILD", str(Path(tempfile.gettempdir())/"spectral-build"))) / 'stata'
def command(*args): return subprocess.check_output(args,text=True)
def inspect(path, minimum):
 text=command('otool','-l',str(path)); records=[]; cmd=None
 for line in text.splitlines():
  fields=line.split()
  if len(fields)==2 and fields[0]=='cmd': cmd=fields[1]
  if len(fields)==2 and ((cmd=='LC_BUILD_VERSION' and fields[0]=='minos') or (cmd=='LC_VERSION_MIN_MACOSX' and fields[0]=='version')):
   assert fields[1]==minimum,(path,fields,minimum)
   records.append(fields[1])
 count=1 if path.suffix!='.a' else sum(not n.startswith('__.SYMDEF') for n in command('ar','-t',str(path)).splitlines())
 assert len(records)==count,(path,len(records),count)
 return count
receipt={}
for plat,arch,minimum in [('macarm64','arm64','11.0'),('macintel64','x86_64','10.15')]:
 for pkg in ['specbunch','specrd']:
  binary=B/'net'/pkg/(pkg+'_'+plat+'.plugin'); tree=B/pkg/plat/'core'
  exports=[s.split()[-1] for s in command('nm','-gU',str(binary)).splitlines() if s.strip()]
  assert sorted(exports)==['_pginit','_stata_call'],exports
  raw=binary.read_bytes(); assert b'__KMP_REGISTERED_LIB_' not in raw
  deps=command('otool','-L',str(binary))
  assert all(s.strip().startswith(('/usr/lib/','/System/Library/')) for s in deps.splitlines()[2:]),deps
  assert command('lipo','-archs',str(binary)).strip()==arch
  inspect(binary,minimum)
  objs=list((tree/'CMakeFiles'/(pkg+'.dir')).rglob('*.o'))
  # Ignore obsolete stplugin.c objects retained by an incremental CMake build.
  objs=[p for p in objs if p.name!='stplugin.c.o']
  assert len(objs)==(8 if pkg=='specbunch' else 4),(pkg,len(objs))
  archives=list(tree.glob('private-openmp-*/lib/libomp_private.a'))
  if pkg=='specbunch': archives+=list(tree.glob('highs/src/highs_ext-build/lib/libhighs.a'))
  count=sum(inspect(p,minimum) for p in objs+archives)
  receipt[pkg+'_'+plat]={'sha256':hashlib.sha256(raw).hexdigest(),'objects':count,'minimum':minimum,'exports':exports,'dependencies':deps,'registry_marker_absent':True}
(B/'audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
