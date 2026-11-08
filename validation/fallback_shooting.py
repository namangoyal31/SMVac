"""Shooting fallback for scan points where the primary solver (C++ Chebyshev collocation) failed.

Failures occur where the bounce crosses the frozen-coupling kink at M_Pl: there the potential is not smooth, and
spectral collocation loses exponential convergence (docs/VALIDATION_REPORT.md section 6). The two independent shooting
solvers M2 and M1 of the Python reference handle the kink through adaptive step-size control. A fallback result is
ACCEPTED only if both converge, agree to < 1e-8 relative, and satisfy the Derrick identity to < 1e-8. The dilatation
identity is reported but degraded by quadrature across the kink.
If both solvers find that EVERY shot undershoots and the Coleman-Glaser-Martin optimum escapes to arbitrarily large
field, the point is classified 'no-finite-bounce': the lambda minimum lies just below M_Pl and lambda barely rises
before the frozen plateau, so the infimum 8 pi^2/(3|lambda(M_Pl)|) is approached only by bubbles escaping into the
plateau (an artefact of freezing the couplings; Planck-sensitive by construction).

Usage: python validation/fallback_shooting.py results/controlled_cw1.csv[.gz] [--workers 10] [--prescription cw1]
Output: <input stem>_fallback.csv next to the input.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))

H0_GEV = 67.4 / 3.0856775814913673e19 * 6.582119569e-25    # Buttazzo et al. eq. (82) convention, as in src/decay.cpp
LN_M_PLANCK = math.log(1.22e19)


def log10_p(S, lnLambdaB):
    return (4.0 * (lnLambdaB - math.log(H0_GEV)) + math.log(0.15) - S) / math.log(10.0)


def solve_point(args):
    Mh, Mt, presc, matching = args
    from smvac.reference import analytic as A, bounce as B, potential as P, sm
    pot = P.sm_potential(sm.Inputs(Mh, Mt, matching=matching), presc)
    an = A.estimates(pot)
    iv = B.admissible_interval(pot, an.x_min, half=12.0)
    row = dict(Mh=Mh, Mt=Mt, S_FL=an.S_FL, S_CGM=an.S_CGM)
    if iv is None:
        return dict(row, fallback_status="no-admissible-interval")
    lo, hi = iv[0], min(iv[1], an.x_min + 8.0)
    r2 = B.solve_m2(pot, lo, hi, n_scan=25)
    r1 = B.solve_m1(pot, lo, hi, n_scan=25)
    row.update(S_M2=r2.S, S_M1=r1.S, derrick_M2=r2.derrick, dilatation_M2=r2.dilatation, derrick_M1=r1.derrick,
               log10_phi0=r2.lnphi0 / math.log(10.0) if r2.ok else float("nan"),
               log10_LambdaB=-r2.lnR / math.log(10.0) if r2.ok else float("nan"))
    agree = r2.ok and r1.ok and abs(r2.S - r1.S) < 1e-8 * abs(r2.S)
    # plateau diagnostics: how far below M_Pl the minimum lies and how much lambda rises before the frozen plateau
    row.update(scan_M2=r2.info.get("scan", ""), scan_M1=r1.info.get("scan", ""), CGM_interior=int(an.CGM_interior),
               log10_mu_min=an.x_min / math.log(10.0), plateau_rise=pot.lam_hi - an.lambda_min,
               S_plateau_infimum=8.0 * math.pi ** 2 / (3.0 * abs(pot.lam_hi)) if pot.lam_hi < 0 else float("inf"))
    if agree and abs(r2.derrick) < 1e-8 and abs(r1.derrick) < 1e-8:
        p = log10_p(r2.S, -r2.lnR)
        row.update(fallback_status="metastable" if p < 0 else "unstable", S_fallback=r2.S,
                   dS_fallback=abs(r2.S - r1.S) + abs(r2.derrick) * abs(r2.S), log10_p=p,
                   planck_sensitive=int(r2.lnphi0 > LN_M_PLANCK))
    elif (not r2.ok and not r1.ok and set(row["scan_M2"]) <= {"u", "s"} and set(row["scan_M1"]) <= {"u", "s"}
          and not an.CGM_interior):
        # every shot undershoots and the best FL trial escapes to arbitrarily large field: no finite bounce. The
        # infimum of the action is approached by bubbles escaping into the frozen-coupling plateau above M_Pl.
        row.update(fallback_status="no-finite-bounce", planck_sensitive=1)
    else:
        row.update(fallback_status="failed")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--prescription", default="cw1")
    ap.add_argument("--matching", default="buttazzo2013")
    a = ap.parse_args()
    df = pd.read_csv(a.scan, comment="#")
    todo = df[df.status == "solver-failure"]
    jobs = [(float(r.Mh), float(r.Mt), a.prescription, a.matching) for r in todo.itertuples()]
    print(f"{len(jobs)} failed points -> shooting fallback with {a.workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        rows = list(ex.map(solve_point, jobs, chunksize=1))
    out = pd.DataFrame(rows)
    stem = a.scan[:-3] if a.scan.endswith(".gz") else a.scan
    path = stem.replace(".csv", "_fallback.csv")
    with open(path, "w", newline="\n") as fh:
        fh.write(f"# shooting fallback (validation/fallback_shooting.py) for solver failures in {os.path.basename(a.scan)}\n")
        fh.write("# acceptance: M2 and M1 converge, |S_M2 - S_M1| < 1e-8 S, |Derrick| < 1e-8 for both\n")
        out.to_csv(fh, index=False)
    print(out.fallback_status.value_counts().to_dict())
    print("wrote", path)


if __name__ == "__main__":
    main()
