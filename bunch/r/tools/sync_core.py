import os, tempfile
from pathlib import Path
import shutil, hashlib, re
p = Path(__file__).resolve().parents[1]
packages = p.parents[1]
src = p / "src"
for name, origin in [("core", p.parent / "core"), ("common", packages / "common")]:
    target = src / name
    if target.exists(): shutil.rmtree(target)
    if name == "core":
        for sub in ("src", "include", "vendor"):
            shutil.copytree(origin/sub, target/sub)
    else:
        target.mkdir()
        for f in ("sb_common.h", "sbc_internal.hpp", "sbc_threads.hpp"): shutil.copy2(origin/f,target/f)
        shutil.copytree(origin/"shim", target/"shim")
# The checked-in trimmed HiGHS tree and translation-unit list are authoritative.
# Refreshing HiGHS is a deliberate dependency update, not needed to build a clone.
rows = []
for root in (src/"core", src/"common", src/"highs"):
    if root.exists():
        for f in sorted(root.rglob("*")):
            if f.is_file(): rows.append(hashlib.sha256(f.read_bytes()).hexdigest()+"  "+f.relative_to(p).as_posix())
(p/"VENDOR.md").write_text("# Byte-copy source SHA-256 manifest\n\nCore/common from sibling directories; HiGHS 1.8.0 from the pinned release.\n\n```\n"+"\n".join(rows)+"\n```\n")

shutil.copy2(p/"VENDOR.md", p/"inst"/"VENDOR.md")

# Keep the installed user reference aligned with the release documentation.
shutil.copy2(packages/"docs"/"REFERENCE.md", p/"inst"/"REFERENCE.md")
