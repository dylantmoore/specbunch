"""Verify vendored SHA-256; optionally compare original research-relative sources."""
import argparse, hashlib, re
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--research-root',type=Path)
    a=p.parse_args();core=Path(__file__).resolve().parents[1]
    rows=re.findall(r'\|\s*`([^`]+)`\s*\|\s*`([^`]+)`\s*\|\s*`([0-9a-f]{64})`\s*\|',(core/'vendor/VENDOR.md').read_text())
    if not rows:raise SystemExit('No hashes in vendor manifest')
    bad=[]
    for filename,origin,digest in rows:
        paths=[core/'vendor'/filename]
        if a.research_root:paths.append(a.research_root/origin)
        ok=all(f.is_file() and hashlib.sha256(f.read_bytes()).hexdigest()==digest for f in paths)
        print('PASS' if ok else 'FAIL',filename)
        if not ok:bad.append(filename)
    print(f'{len(rows)} files; {len(bad)} failures; mode: '+('SHA + research' if a.research_root else 'SHA only'))
    return bool(bad)
if __name__=='__main__':raise SystemExit(main())
