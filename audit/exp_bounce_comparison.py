"""E6: legacy S_exact / S_approx vs independent rigorous bounds and exact shooting bounce,
for the SAME legacy physics inputs (matching, betas, lambda_eff with mu = phi).
Two potentials: pure SM lambda_eff (what S_approx uses) and lambda_eff + phi^6/(6 M_Pl^2)
(what the legacy S_exact integrand uses).  Output: audit/results/bounce_comparison.csv"""
import math, csv, sys
import numpy as np
from audit import legacy_port as L, independent as I

POINTS = [(125.1, 173.1), (120.0, 172.0), (130.0, 178.0), (115.0, 180.0), (134.75, 176.5)]
if len(sys.argv) > 1:
    POINTS = [tuple(map(float, p.split(','))) for p in sys.argv[1:]]
rows = []
for Mh, Mt in POINTS:
    leg = L.classify(Mh, Mt)
    rg = I.RGSolution(Mh, Mt)
    for phi6 in (False, True):
        fpot = I.FieldPotential(rg, include_phi6=phi6)
        b = I.tree_bounds(fpot)
        row = dict(Mh=Mh, Mt=Mt, phi6=phi6, legacy_status=leg['status'], legacy_S_exact=leg['S_exact'],
                   legacy_S_approx=leg['S_approx'], lambda_min=b['lambda_min'],
                   log10_phi_at_min=b['lnphi_min'] / math.log(10), S_lower=b['S_lower'], S_upper=b['S_upper'],
                   S_shoot=float('nan'), derrick=float('nan'), dilatation=float('nan'), log10_phi0=float('nan'), note='')
        c = b['lnphi_min']
        tp = I.TabulatedPotential(fpot, min(fpot.t_lo / 2, c - 30), c + 30)
        br, kinds = I.find_shooting_bracket(tp, c - 6, c + 6, n=25)
        if br:
            s = I.shoot_bounce(tp, *br[0])
            row.update(S_shoot=s['S'], derrick=s['derrick'], dilatation=s['dilatation'],
                       log10_phi0=s['lnphi0'] / math.log(10), note=f'{len(br)} bracket(s)')
        else:
            row['note'] = 'no under/over bracket: ' + ''.join(k[1][0] for k in kinds)
        rows.append(row)
        print({k: (f"{v:.7g}" if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
out = 'audit/results/bounce_comparison.csv' if len(sys.argv) == 1 else 'audit/results/bounce_comparison_extra.csv'
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('wrote', out)
