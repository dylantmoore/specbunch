#!/usr/bin/env python3
import os, tempfile
"""Build outside Drive, stage one net-install directory; never install globally."""
import argparse, os, subprocess, shutil
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--platform',choices=['macarm64','macintel64','linux64','win64'],required=True)
p.add_argument('--build-root',type=Path,default=Path(os.environ.get("SPECPKG_BUILD", str(Path(tempfile.gettempdir())/"spectral-build"))) / 'stata')
a=p.parse_args(); src=Path(__file__).resolve().parent; pkg='specbunch'
b=a.build_root/pkg/a.platform; b.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,TMPDIR=str(a.build_root/'tmp'))
(a.build_root/'tmp').mkdir(parents=True,exist_ok=True)
args=['cmake','-S',str(src),'-B',str(b),'-DCMAKE_BUILD_TYPE=Release','-DSTATA_PLATFORM='+a.platform]
if a.platform.startswith('mac'):
 arch='arm64' if a.platform=='macarm64' else 'x86_64'
 args+=['-DCMAKE_OSX_ARCHITECTURES='+arch,'-DCMAKE_OSX_DEPLOYMENT_TARGET='+('11.0' if arch=='arm64' else '10.15'),'-DSBC_PRIVATE_OPENMP=ON']
subprocess.run(args,check=True,env=env)
subprocess.run(['cmake','--build',str(b),'--target','stata_plugin','--parallel','2',*(['--verbose'] if os.name=='nt' else [])],check=True,env=env)
stage=a.build_root/'net'/pkg; stage.mkdir(parents=True,exist_ok=True)
for f in src.iterdir():
 if f.suffix in ['.ado','.sthlp','.pkg','.toc','.txt']: (stage/f.name).write_bytes(f.read_bytes())
shutil.copy2(b/'core'/(pkg+'_'+a.platform+'.plugin'),stage)
# A platform-limited local manifest must not reference unbuilt platforms.
lines=(src/(pkg+'.pkg')).read_text().splitlines()
(stage/(pkg+'.pkg')).write_text('\n'.join(x for x in lines if not x.endswith('.plugin') or (stage/x[2:]).exists())+'\n')

subprocess.run([os.environ.get("PYTHON", "python"), "-B", str(src.parents[1]/"tools/source_archive.py"),
                "--repo", str(src.parents[1]), "--output", str(stage/"specbunch_source.tar.gz"),
                "--cache", str(a.build_root/"source-cache"), "--stata-build", str(b)], check=True)
