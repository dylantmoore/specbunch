"""Compare encoded outputs, enforcing exact discrete decisions on every platform."""
import math, re

def compare(actual, expected, *, exact=False):
    assert actual['status'] == expected['status'], ('status', actual['status'], expected['status'])
    assert actual['strings'] == expected['strings'], ('strings/decisions', actual['strings'], expected['strings'])
    count = bitwise = 1 + len(expected['strings'])
    maxrel = 0.0
    for group in ('scalars','vectors'):
        assert actual[group].keys() == expected[group].keys(), (group, 'keys')
        for key, want in expected[group].items():
            got = actual[group][key]
            aa, bb = ([got],[want]) if group == 'scalars' else (got,want)
            assert len(aa)==len(bb), (key,'length')
            discrete = bool(re.search(r'status|decision|flag|tier|mask|accept|refus|selected|index|^n_|^n$|weighted|rule|answered|fallback|_ok$|policy_mode|(?:^|[._])rank$', key, re.I))
            for a,b in zip(aa,bb):
                count += 1
                if a['bits']==b['bits']:
                    bitwise += 1
                    continue
                x,y = float.fromhex(a['hex']),float.fromhex(b['hex'])
                if math.isnan(x) and math.isnan(y) and not exact: continue
                rel = abs(x-y)/max(abs(x),abs(y),1e-300) if math.isfinite(x) and math.isfinite(y) else math.inf
                maxrel=max(maxrel,rel)
                assert not exact and not discrete and rel <= 1e-10, (key,x,y,rel,'exact' if exact or discrete else 'tolerance')
    return dict(compared=count, bitwise=bitwise, maxrel=maxrel)
