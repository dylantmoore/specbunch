"""Create corresponding source/relinking materials for the LGPL component.
The archive is generated outside the checkout and shipped with binary artifacts.
"""
import argparse, hashlib, io, shutil, tarfile, urllib.request
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--cache',type=Path,required=True);p.add_argument('--stata-build',type=Path);a=p.parse_args()
r=a.repo.resolve();a.cache.mkdir(parents=True,exist_ok=True);a.output.parent.mkdir(parents=True,exist_ok=True)
deps=[('HiGHS-1.8.0.tar.gz','https://github.com/ERGO-Code/HiGHS/archive/refs/tags/v1.8.0.tar.gz','e184e63101cf19688a02102f58447acc7c021d77eef0d3475ceaceb61f035539'),('openmp-21.1.8.src.tar.xz','https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.8/openmp-21.1.8.src.tar.xz','856b023748b41ac7b2c83fd8e9f765ff48a4df2fe6777d2811ef7c7ed8f2f977'),('cmake-21.1.8.src.tar.xz','https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.8/cmake-21.1.8.src.tar.xz','85735f20fd8c81ecb0a09abb0c267018475420e93b65050cc5b7634eab744de9')]
for name,url,digest in deps:
 f=a.cache/name
 if not f.exists():
  temp=f.with_suffix(f.suffix+'.tmp')
  with urllib.request.urlopen(url,timeout=120) as response,temp.open('wb') as out:shutil.copyfileobj(response,out)
  temp.replace(f)
 if hashlib.sha256(f.read_bytes()).hexdigest()!=digest:raise SystemExit('Source archive hash mismatch: '+name)
with tarfile.open(a.output,'w:gz') as out:
 for rel in ['LICENSE','NOTICE','THIRD_PARTY_NOTICES.md','licenses','docs/RELINKING.md','bunch/core','common','bunch/stata','bunch/python','bunch/julia','tools/source_archive.py']:
  src=r/rel
  files=sorted(src.rglob('*')) if src.is_dir() else [src]
  for file in files:
   if not file.is_file() or any(x in file.relative_to(r).parts for x in ['__pycache__','.pytest_cache','build','.git']):continue
   # Explicit reads hydrate Drive placeholders before their size is recorded.
   data=file.read_bytes()
   info=out.gettarinfo(str(file),arcname='specbunch-source/'+file.relative_to(r).as_posix())
   info.size=len(data)
   out.addfile(info,io.BytesIO(data))
 if a.stata_build:
  for stem in ['plugin.cpp','entry.c']:
   candidates=[p for p in a.stata_build.rglob('*') if p.name in [stem+'.o',stem+'.obj'] and 'specbunch.dir' in p.parts]
   if len(candidates)!=1:raise SystemExit('Expected one Stata wrapper object: '+str(candidates))
   out.add(candidates[0],arcname='specbunch-source/prebuilt-stata/'+candidates[0].name)
 for name,_,_ in deps:out.add(a.cache/name,arcname='specbunch-source/third_party/'+name)
print(a.output, a.output.stat().st_size)
