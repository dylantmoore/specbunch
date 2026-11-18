import os, tempfile
"""Check shipped Mach-O minos AND each linked core/static-archive object.

A dylib's LC_BUILD_VERSION alone cannot certify its input objects. This audit
also inspects the exact core .o files and HiGHS/private-libomp archive members.
"""
import argparse, hashlib, json, re, subprocess, tempfile, zipfile
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--root',type=Path,default=Path(os.environ.get("SPECPKG_BUILD", str(Path(tempfile.gettempdir())/"spectral-build"))))
p.add_argument('--arch',choices=['arm64','x86_64'],default='arm64')
p.add_argument('--build-root',type=Path)
p.add_argument('--wheel-dir',type=Path)
a=p.parse_args()
build=a.build_root or a.root/'py/build'; wheels=a.wheel_dir or a.root/'py/dist'
minimum='11.0' if a.arch=='arm64' else '10.15'
tag='macosx_'+minimum.replace('.','_')+'_'+a.arch

def inspect(path):
    text=subprocess.check_output(['otool','-l',str(path)],text=True)
    member=str(path); command=None; records=[]
    for line in text.splitlines():
        if line.endswith(':') and not line.startswith((' ','Archive :')):
            member=line[:-1]
        fields=line.split()
        if len(fields)==2 and fields[0]=='cmd': command=fields[1]
        if len(fields)==2 and ((command=='LC_BUILD_VERSION' and fields[0]=='minos') or
                              (command=='LC_VERSION_MIN_MACOSX' and fields[0]=='version')):
            assert fields[1]==minimum,(member,fields[1],minimum)
            records.append({'object':member,'command':command,'minos':fields[1]})
    count=1
    if path.suffix=='.a':
        members=subprocess.check_output(['ar','-t',str(path)],text=True).splitlines()
        count=sum(not n.startswith('__.SYMDEF') for n in members)
    assert len(records)==count,(path,len(records),count,'missing/duplicate minimum-version command')
    return records

result={'architecture':a.arch,'minimum_macos':minimum,'packages':{}}
for name in ['specbunch','specrd']:
    candidates=list(wheels.glob(name+'-*-py3-none-'+tag+'.whl'))
    assert len(candidates)==1,candidates
    wheel=candidates[0]; package={'wheel':wheel.name,'wheel_sha256':hashlib.sha256(wheel.read_bytes()).hexdigest()}
    with zipfile.ZipFile(wheel) as z, tempfile.TemporaryDirectory(dir=a.root/'py/tmp') as temp:
        metadata=z.read(next(n for n in z.namelist() if n.endswith('.dist-info/WHEEL'))).decode()
        assert 'Tag: py3-none-'+tag in metadata,metadata
        library=next(n for n in z.namelist() if n.endswith('.dylib'))
        path=Path(temp)/Path(library).name; path.write_bytes(z.read(library))
        assert subprocess.check_output(['lipo','-archs',str(path)],text=True).strip()==a.arch
        package['shipped_dylib']=inspect(path)
        package['library_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    tree=build/name/'core'
    core_objects=sorted((tree/'CMakeFiles'/(name+'.dir')).rglob('*.o'))
    assert len(core_objects)==(6 if name=='specbunch' else 2),core_objects
    package['core_objects']=[r for path in core_objects for r in inspect(path)]
    archives=sorted(tree.glob('private-openmp-*/lib/libomp_private.a'))
    if name=='specbunch': archives+=list(tree.glob('highs/src/highs_ext-build/lib/libhighs.a'))
    assert len(archives)==(2 if name=='specbunch' else 1),archives
    package['archive_objects']=[r for path in archives for r in inspect(path)]
    package['object_count']=len(package['core_objects'])+len(package['archive_objects'])
    result['packages'][name]=package
    print(name,tag,'shipped dylib and',package['object_count'],'objects PASS',flush=True)
out=a.root/('py/deployment/audit_'+a.arch+'.json')
out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,indent=2)+'\n')
