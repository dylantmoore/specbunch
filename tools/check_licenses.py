"""Check licence propagation and built distribution contents (stdlib only)."""
from pathlib import Path
import argparse, re, tarfile, zipfile
p=argparse.ArgumentParser();p.add_argument('--wheel',type=Path);p.add_argument('--r-tarball',type=Path);p.add_argument('--stata',type=Path);a=p.parse_args()
def legal(data):return data.replace(b'\r\n',b'\n')
r=Path(__file__).resolve().parents[1];name='specbunch' if (r/'bunch').is_dir() else 'specrd';f='bunch' if name=='specbunch' else 'rd';j='SpecBunch' if f=='bunch' else 'SpecRD'
for dest in [r/f/'python',r/f/'r/inst',r/f/'julia'/j]:
 for fn in ['LICENSE','NOTICE','THIRD_PARTY_NOTICES.md','CITATION.cff']:
  assert legal((dest/fn).read_bytes())==legal((r/fn).read_bytes()),str(dest/fn)
 for source in (r/'licenses').iterdir():assert legal((dest/'licenses'/source.name).read_bytes())==legal(source.read_bytes()),str(dest/source.name)
assert (r/f/'r/LICENSE').read_text()=='YEAR: 2026\nCOPYRIGHT HOLDER: Dylan Moore\n'
assert 'License: MIT + file LICENSE' in (r/f/'r/DESCRIPTION').read_text()
assert not list(r.rglob('stplugin.c')) and not list(r.rglob('stplugin.h')), 'SDK source must not be redistributed'
for fn in ['libomp-private.patch','libomp-private-windows.patch']:assert (r/'licenses'/fn).read_bytes()==(r/'common/third_party'/fn).read_bytes()
if a.wheel:
 with zipfile.ZipFile(a.wheel) as z:
  names=z.namelist();prefix=name+'-0.1.0.dist-info/licenses/'
  for fn in ['LICENSE','NOTICE','THIRD_PARTY_NOTICES.md']+['licenses/'+p.name for p in (r/'licenses').iterdir()]:
   assert prefix+fn in names,fn
   assert legal(z.read(prefix+fn))==legal((r/fn).read_bytes()), 'stale wheel notice: '+fn
  if name=='specbunch':assert name+'/licenses/corresponding-source.tar.gz' in names
  if 'manylinux' in a.wheel.name:
   info=z.read(prefix+'licenses/GCC-CORRESPONDING-SOURCE.txt').decode()
   digest=re.search(r'^SHA-256: ([0-9a-f]{64})$',info,re.M).group(1)
   bundled=[n for n in names if re.search(r'libgomp-[0-9a-f]{8}\.so',n)]
   assert len(bundled)==1 and ('libgomp-'+digest[:8]+'.so') in bundled[0], 'GCC source notice does not identify bundled runtime'
if a.r_tarball:
 with tarfile.open(a.r_tarball) as t:
  names=t.getnames()
  assert legal(t.extractfile(name+'/inst/COPYRIGHTS').read())==legal((r/'THIRD_PARTY_NOTICES.md').read_bytes()), 'stale R notice'
  for fn in ['LICENSE','inst/COPYRIGHTS','inst/CITATION','inst/NOTICE']+['inst/licenses/'+p.name for p in (r/'licenses').iterdir()]:assert name+'/'+fn in names,fn
  for source in (r/'licenses').iterdir():assert legal(t.extractfile(name+'/inst/licenses/'+source.name).read())==legal(source.read_bytes()), 'stale R licence: '+source.name
if a.stata:
 assert legal((a.stata/(name+'_THIRD_PARTY_NOTICES.txt')).read_bytes())==legal((r/'THIRD_PARTY_NOTICES.md').read_bytes())
 for row in (a.stata/(name+'.pkg')).read_text().splitlines():
  if row.startswith('f '):assert (a.stata/row[2:]).is_file(),row
 for fn in ['LICENSE','NOTICE','THIRD_PARTY_NOTICES','CITATION']:assert (a.stata/(name+'_'+fn+'.txt')).is_file()
 for source in (r/'licenses').iterdir():assert legal((a.stata/(name+'_'+source.stem+'.txt')).read_bytes())==legal(source.read_bytes()), 'stale Stata licence: '+source.name
 if name=='specbunch':assert (a.stata/'specbunch_source.tar.gz').is_file()
print('PASS:',name,'licences, notices, citations and requested artifact contents')
