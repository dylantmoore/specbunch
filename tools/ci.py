"""Small cross-platform CI entry points. All persistent outputs stay in --root."""
import argparse, os, subprocess, sys, shutil, platform
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('job',choices=['cores','fixtures','stata','julia','wheels','wheel-test'])
p.add_argument('--root',type=Path,required=True);p.add_argument('--platform');p.add_argument('--family',choices=['bunch']);a=p.parse_args()
r=Path(__file__).resolve().parents[1];b=a.root.resolve();b.mkdir(parents=True,exist_ok=True)
os.environ.update(SPECPKG_BUILD=str(b),PYTHONDONTWRITEBYTECODE='1',CMAKE_BUILD_PARALLEL_LEVEL='3')
if platform.system()=='Windows':os.environ.update(CC='clang-cl',CXX='clang-cl',CMAKE_GENERATOR='Ninja')
def run(*args,**kw):subprocess.run(list(map(str,args)),check=True,**kw)
def lib(name):return ('' if os.name=='nt' else 'lib')+name+('.dll' if os.name=='nt' else '.dylib' if sys.platform=='darwin' else '.so')
def build(src,dest,*opts):
    run('cmake','-S',r/src,'-B',dest,'-DCMAKE_BUILD_TYPE=Release',*opts)
    run('cmake','--build',dest,'--config','Release','--parallel','3',*(['--verbose'] if os.name=='nt' else []))
run(sys.executable, '-B', r/'tools/check_common.py')
families=[a.family] if a.family else ['bunch']
if a.job in ('cores','fixtures'):
    for f in families:build(f+'/core',b/'cores'/f)
    if a.job=='fixtures':
        run('cmake','--build',b/'cores/bunch','--config','Release','--target','sb_test_joint_voigt')
        run(b/'cores/bunch'/('sb_test_joint_voigt.exe' if os.name=='nt' else 'sb_test_joint_voigt'))
        run(sys.executable,'-B',r/'tests/xlang/generate.py','--bunch',b/'cores/bunch'/lib('specbunch'),'--root',b/'xlang/fixtures')
elif a.job=='wheels':
    for f in families:
        name='specbunch' if f=='bunch' else 'specrd'
        run(sys.executable,'-B',r/f/'python/tools/build.py','--stage-only','--root',b)
        run(sys.executable,'-m','cibuildwheel','.', '--output-dir',b/'py/dist',cwd=b/'py/stage'/name)
        for wheel in (b/'py/dist').glob('*.whl'):run(sys.executable,'-B',r/'tools/check_licenses.py','--wheel',wheel)
elif a.job=='wheel-test':
    run(sys.executable,'-m','pip','install',*sorted((b/'py/dist').glob('*.whl')),'pytest')
    (b/'py').mkdir(exist_ok=True)
    run(sys.executable,'-B',r/'tests/xlang/run_python.py','--root',b,'--worker','both','--tolerance')
    # Linux produces the shared oracle. Each native wheel must also match its own C ABI bit-for-bit.
    if sys.platform=='darwin' and platform.machine()=='arm64':
        import specbunch
        local=b/'native-oracle';(local/'py').mkdir(parents=True,exist_ok=True)
        def locate(mod,name):return next(Path(mod.__file__).parent.rglob(lib(name)))
        run(sys.executable,'-B',r/'tests/xlang/generate.py','--bunch',locate(specbunch,'specbunch'),'--root',local/'xlang/fixtures')
        run(sys.executable,'-B',r/'tests/xlang/run_python.py','--root',local,'--worker','both')
elif a.job=='stata':
    build('tests/no_stata',b/'host')
    for f in families:
        name='specbunch' if f=='bunch' else 'specrd'
        run(sys.executable,'-B',r/f/'stata/build.py','--platform',a.platform,'--build-root',b/'stata')
        run(sys.executable,'-B',r/'tools/check_licenses.py','--stata',b/'stata/net'/name)
        run(sys.executable,'-B',r/'tests/no_stata/run.py','--host',b/'host'/lib('stata_host'),
            '--plugin',b/'stata'/name/a.platform/'core'/(name+'_'+a.platform+'.plugin'),
            '--fixtures',b/'xlang/fixtures','--family',f,'--output',b/(name+'-plugin.json'))
elif a.job=='julia':
    j=b/'julia';os.environ.update(SB_JULIA_BUILD=str(j),JULIA_DEPOT_PATH=str(j/'depot'))
    for f in families:build(f+'/core',j/'build'/f)
    run('julia','--startup-file=no',r/'bunch/julia/tools/install_local.jl',*families)
