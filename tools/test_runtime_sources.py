"""Regression: a GCC linker script must not be mistaken for bundled libgomp."""
from pathlib import Path
import hashlib, os, runpy, shutil, tempfile
from unittest.mock import patch
from types import SimpleNamespace
repo=Path(__file__).resolve().parents[1]
family='bunch' if (repo/'bunch').exists() else 'rd'
with tempfile.TemporaryDirectory() as d:
 root=Path(d).resolve();(root/'tools').mkdir();(root/'licenses').mkdir()
 script=root/'tools/runtime_sources.py';shutil.copyfile(repo/family/'python/tools/runtime_sources.py',script)
 elf=root/'libgomp.so.1.0.0';elf.write_bytes(b'\x7fELFactual-runtime')
 linker=root/'libgomp.so';linker.write_text('INPUT ( libgomp.so.1 -lgomp_nonshared )')
 original_read=Path.read_text
 def read(path,*args,**kw):
  if str(path)=='/proc/self/maps':return '0000-ffff r-xp 0000 00:00 1 '+str(elf)+'\n'
  return original_read(path,*args,**kw)
 def command(args,**kw):
  if args[0]=='gcc':return str(linker if args[1].endswith('=libgomp.so') else root/'missing.a')+'\n'
  source='gcc-8.5.0.src.rpm' if args[-1]==str(elf) else 'gcc-toolset-14-gcc.src.rpm'
  return source if '--qf' in args else source.removesuffix('.src.rpm')+'.x86_64\n'
 class Response:
  status=200
  def __enter__(self):return self
  def __exit__(self,*args):pass
 with patch.dict(os.environ,{'CC':'gcc'}),patch('ctypes.CDLL',return_value=object()),patch.object(Path,'read_text',read),patch('subprocess.check_output',command),patch('urllib.request.urlopen',return_value=Response()):
  runpy.run_path(str(script),run_name='__main__')
 text=(root/'licenses/GCC-CORRESPONDING-SOURCE.txt').read_text()
 assert str(elf) in text.split('Compiler link input')[0]
 assert str(linker) in text.split('Compiler link input')[1]
 assert 'SHA-256: '+hashlib.sha256(elf.read_bytes()).hexdigest() in text.split('Compiler link input')[0]
 assert 'gcc-8.5.0.src.rpm' in text.split('Compiler link input')[0]
 assert 'gcc-toolset-14-gcc.src.rpm' in text.split('Compiler link input')[1]
print('PASS: actual ELF runtime and compiler linker script receive distinct source provenance')
