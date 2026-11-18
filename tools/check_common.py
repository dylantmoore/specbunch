"""Verify the recorded SHA-256 of every common source file."""
from pathlib import Path
import hashlib,json,sys
root=Path(__file__).resolve().parents[1]
actual={str(p.relative_to(root/'common')):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((root/'common').rglob('*')) if p.is_file()}
expected=json.loads((root/'tools/common.sha256.json').read_text())
if actual != expected:
    raise SystemExit('common source drift: update manifest deliberately and compare both release repositories')
print('PASS: common source hashes',len(actual))
