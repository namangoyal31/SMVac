"""Check that audit/legacy_port.py reproduces the compiled legacy C++ solver.

Usage: python -m audit.validate_port <legacy_driver_csv>
where the CSV was produced by the legacy driver (audit/legacy_driver.cpp) for the same points.
"""
import sys, csv, math
from audit import legacy_port as L

rows = list(csv.DictReader(open(sys.argv[1])))
worst = 0.0
for r in rows:
    Mh, Mt = float(r['Mh']), float(r['Mt'])
    res = L.classify(Mh, Mt)
    ok = res['status'] == int(r['Stability'])
    for k in ('S_exact', 'S_approx'):
        a, b = res[k], float(r[k])
        rel = abs(a - b) / max(abs(b), 1e-300)
        worst = max(worst, rel)
        ok = ok and rel < 1e-9
    print(f"Mh={Mh:8.3f} Mt={Mt:8.3f}  C++: st={r['Stability']} S_exact={float(r['S_exact']):.12g} S_approx={float(r['S_approx']):.12g} | port: st={res['status']} S_exact={res['S_exact']:.12g} S_approx={res['S_approx']:.12g} {'OK' if ok else 'MISMATCH'}")
print("max relative difference:", worst)
