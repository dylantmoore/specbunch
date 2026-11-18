# Vendored source provenance

These unmodified bytes are verified against the SHA-256 table below. Original
research-relative paths identify provenance only and are not build dependencies.
Run `tools/check_vendor.py` for SHA-only validation; optional
`--research-root PATH` also compares the original sources.

| vendored file | original research path | sha256 |
|---|---|---|
| `llcore.hpp` | `cpp/part1/src/llcore.hpp` | `07f124e5ec230d887361c13ddff9741201c4d5613be6bd592a5d2ac8e5572590` |
| `p1lad.hpp` | `cpp/part1/src/puf/p1lad.hpp` | `336571293292d2d852498d5b4f7b704eed6da1ab1f41aa128befdd1006936d45` |
| `p1w.hpp` | `cpp/part1/src/puf/p1w.hpp` | `a1034929cd8157205303361afbf88c25954ca31cfc68bb04306def053fd6c736` |
| `c2/kernel_puf.cpp` | `cpp/part1/src/puf/c2/kernel_puf.cpp` | `b5721dc2b9e90f287d64b75d5506aba73623f2eb9e2b09dcc9183864adce907d` |
| `c2/probit.cpp` | `cpp/part1/src/puf/c2/probit.cpp` | `317bf165afd8f2f4e4d3aa1175cf9ef5653525548ec2949b110719ef2db535de` |
| `famfit.hpp` | `cpp/part1/src/famfit.hpp` | `cb810bbff9264d56005fa5f0e6b6a34b36a670b537a24ae80707c56103b91bef` |
| `puf/p1famv.hpp` | `cpp/part1/src/puf/p1famv.hpp` | `be087caa3a53ae2c385c95646160fc237e80c5cfff0d543518d7fb1b3d7d87d8` |
| `cfj/cfj.cpp` | `cpp/part1/cfj/cfj.cpp` | `df7efb49591ad1f660fbe21211a3ff010f0494c988c3eaad4600dee57185cf8d` |
| `p2/p2core.hpp` | `cpp/part2/p2core.hpp` | `fec1439c08ec4086e5fb1d6f59a402bbee4ac8c7838e3845a62a1c4d728b74af` |
| `p2/p2pilot.hpp` | `cpp/part2/p2pilot.hpp` | `e99f06b5f885e226c2f28bc9b50eebb1c38e0ba347ee00d086fb25297ad98891` |

The non-vendored src/mean_engine.hpp adapts p2/p2core.hpp to retain ranked bands and use a coherent fallback. The original vendored source remains unchanged.
