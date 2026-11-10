"""Convergence and sensitivity study of the numerical reference (Stage D evidence).

For representative cases (thick wall, onset of thin wall, thin wall, SM near-critical):
  * M4 Chebyshev collocation: S(N) for N = 32 ... 1024 and the identity residuals;
  * M2 and M1 shooting: S(rtol) for rtol = 1e-6 ... 1e-13;
  * SM cases: sensitivity to the RG integration tolerance and to the Chebyshev degree of lambda(x).
The reference value per case is the M1/M2/M4 consensus at their tightest settings.
Output: validation/results/convergence.csv
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
from smvac.reference import analytic as A, bounce as B, potential as P, sm  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "results", "convergence.csv")


def cases():
    yield "parabola q=250 (thick)", P.parabola_model(-0.0125, 5e-5, 39.0), (25, 55)
    yield "parabola q=2.5 (thin onset)", P.parabola_model(-2.5 * 5e-5, 5e-5, 39.0), (25, 55)
    yield "parabola q=0.1 (thin wall)", P.parabola_model(-0.1 * 5e-5, 5e-5, 39.0), (25, 55)
    yield "SM cw1 (125.15,173.34)", P.sm_potential(sm.Inputs(125.15, 173.34), "cw1"), None
    yield "SM cw1 (125.15,171.149)", P.sm_potential(sm.Inputs(125.15, 171.119296 + 0.03), "cw1"), None


def main():
    rows = []
    for name, pot, search in cases():
        an = A.estimates(pot, search)
        lnh = an.lnh_CGM if math.isfinite(an.lnh_CGM) else an.x_min + 5.0 / 6.0
        iv = B.admissible_interval(pot, an.x_min)
        lo, hi = iv[0], min(iv[1], an.x_min + 8.0)
        for N in (32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024):
            t = time.perf_counter()
            r = B.solve_m4(pot, an.lambda_min, lnh, N=N)
            rows.append(dict(case=name, method="M4", parameter="N", value=N, ok=r.ok, S=r.S, derrick=r.derrick,
                             dilatation=r.dilatation, seconds=time.perf_counter() - t))
            print(f"{name:28s} M4 N={N:5d} ok={r.ok!s:5} S={r.S:.15g} derrick={r.derrick:.1e}", flush=True)
        for m, fn in (("M2", B.solve_m2), ("M1", B.solve_m1)):
            for rt in (1e-6, 1e-8, 1e-10, 1e-11, 1e-12, 1e-13):
                t = time.perf_counter()
                r = fn(pot, lo, hi, n_scan=25, rtol=rt)
                rows.append(dict(case=name, method=m, parameter="rtol", value=rt, ok=r.ok, S=r.S, derrick=r.derrick,
                                 dilatation=r.dilatation, seconds=time.perf_counter() - t))
                print(f"{name:28s} {m} rtol={rt:.0e} ok={r.ok!s:5} S={r.S:.15g} derrick={r.derrick:.1e}", flush=True)
        if name.startswith("SM"):
            Mt = 173.34 if "173.34" in name else 171.119296 + 0.03
            inp = sm.Inputs(125.15, Mt)
            for label, kw in (("rge_rtol=1e-10", dict(rge=1e-10)), ("rge_rtol=1e-13", dict(rge=1e-13)),
                              ("cheb_deg=128", dict(deg=128)), ("cheb_deg=256", dict(deg=256)),
                              ("cheb_deg=512", dict(deg=512))):
                rg = P.RGSolution(inp, rtol=kw.get("rge", 1e-13))
                p2 = P.sm_potential(inp, "cw1", rg=rg, deg=kw.get("deg"))
                a2 = A.estimates(p2)
                r = None
                for N in (256, 512, 1024):
                    r = B.solve_m4(p2, a2.lambda_min, a2.lnh_CGM, N=N)
                    if r.ok and abs(r.derrick) < 1e-9 and abs(r.dilatation) < 1e-9:
                        break
                rows.append(dict(case=name, method="M4", parameter="potential", value=label, ok=r.ok, S=r.S,
                                 derrick=r.derrick, dilatation=r.dilatation, seconds=float("nan")))
                print(f"{name:28s} potential {label:16s} S={r.S:.15g} S_FL={a2.S_FL:.15g}", flush=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
