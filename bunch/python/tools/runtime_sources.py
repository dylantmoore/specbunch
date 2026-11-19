"""Record source access for the ELF runtime actually redistributed by auditwheel.
Run inside manylinux before wheel construction. GCC's libgomp.so can be a
linker script owned by a different source RPM; it is not the loaded runtime.
"""
from pathlib import Path
import ctypes, hashlib, os, subprocess, urllib.request
root=Path(__file__).resolve().parents[1]
def run(*args):return subprocess.check_output(args,text=True).strip()
loaded=ctypes.CDLL('libgomp.so.1')
paths={Path(line.split()[-1]).resolve() for line in Path('/proc/self/maps').read_text().splitlines() if '/libgomp.so.' in line}
if len(paths)!=1:raise SystemExit('Cannot identify one loaded libgomp runtime: '+str(paths))
gomp=paths.pop()
if not gomp.read_bytes().startswith(b'\x7fELF'):raise SystemExit('Runtime is not an ELF library: '+str(gomp))
inputs=[('Redistributed ELF runtime',gomp)]
for filename in ['libgomp.so','libgomp_nonshared.a']:
 path=Path(run(os.environ.get('CC','gcc'),'-print-file-name='+filename)).resolve()
 if path.is_file() and path!=gomp:inputs.append(('Compiler link input (not the shared runtime)',path))
urls={};records=[]
for label,path in inputs:
 source=run('rpm','-qf','--qf','%{SOURCERPM}',str(path))
 if not source.endswith('.src.rpm'):raise SystemExit('No matching GCC source RPM: '+source)
 if source not in urls:
  for prefix in ['vault/8.10','8']:
   for repo in ['BaseOS','AppStream']:
    url=f'https://repo.almalinux.org/{prefix}/{repo}/Source/Packages/{source}'
    try:
     with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=30) as response:
      if response.status==200:urls[source]=url;break
    except Exception:pass
   if source in urls:break
  if source not in urls:raise SystemExit('Matching GCC source unavailable; do not ship: '+source)
 records.append(f'{label}: {path}\nSHA-256: {hashlib.sha256(path.read_bytes()).hexdigest()}\nRPM: {run("rpm","-qf",str(path))}\nSource RPM: {source}\nCorresponding source: {urls[source]}\n')
(root/'licenses/GCC-CORRESPONDING-SOURCE.txt').write_text(
 'The Linux wheel redistributes the ELF runtime below, renamed by auditwheel.\n'
 'GPL-3.0-or-later WITH GCC-exception-3.1.\n\n'+'\n'.join(records)+
 '\nSource RPMs include GCC source, distributor patches and the build specification.\n'
 'Download without charge, extract using rpm2cpio and cpio, and follow the RPM spec.\n'
 'Compiler link inputs are recorded separately to cover any statically used support.\n'
 'Application code remains under its own licence through the GCC runtime exception.\n')
print('Matching runtime source:',records[0])
