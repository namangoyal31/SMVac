"""Systematic comparison of candidate bounce solvers (docs/STAGE_C_PLAN.md section 6).

Runs M0-M4 on the benchmark set B1-B5, then convergence sweeps, and writes
    validation/results/method_comparison.csv      one row per (benchmark, method)
    validation/results/method_convergence.csv     refinement sweeps
Usage: python validation/method_comparison.py [--quick]
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
from scipy.optimize import brentq  # noqa: E402

from smvac.reference import analytic as A, bounce as B, potential as P, sm  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(OUT, exist_ok=True)
QUICK = "--quick" in sys.argv


def sm_case(Mh, Mt, presc="cw1"):
    pot = P.sm_potential(sm.Inputs(Mh, Mt), presc)
    return pot, None


def critical_Mt(Mh, presc="cw1", lo=169.0, hi=174.0):
    """Top mass at which min lambda_P = 0 (absolute-stability boundary) for given Mh."""
    def g(Mt):
        pot = P.sm_potential(sm.Inputs(Mh, Mt), presc)
        return A.global_min(pot)[0]
    return brentq(g, lo, hi, xtol=1e-6)


CSV_PATH = os.path.join(OUT, "method_comparison.csv")
FIELDS = ["case", "method", "ok", "S", "derrick", "dilatation", "lnphi0", "lnR", "seconds", "lambda_min", "x_min",
          "lambda2", "q", "S_FL", "S_curv", "S_CGM", "notes", "info"]


def done_pairs():
    if not os.path.exists(CSV_PATH):
        return set()
    with open(CSV_PATH) as fh:
        return {(r["case"], r["method"]) for r in csv.DictReader(fh)}


def append(row):
    new = not os.path.exists(CSV_PATH)
    with open(CSV_PATH, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


def run_all(name, pot, search, methods=("M0", "M1", "M2", "M3", "M4"), notes=""):
    done = done_pairs()
    methods = [m for m in methods if (name, m) not in done]
    if not methods:
        print(f"  {name}: already done", flush=True)
        return []
    an = A.estimates(pot, search)
    rows = []
    base = dict(case=name, lambda_min=an.lambda_min, x_min=an.x_min, lambda2=an.lambda2,
                q=abs(an.lambda_min) / an.lambda2 if an.lambda2 and an.lambda2 > 0 else float("nan"),
                S_FL=an.S_FL, S_curv=an.S_curv, S_CGM=an.S_CGM, notes=notes + ";".join(an.notes))
    if an.lambda_min >= 0:
        return [dict(base, method="-", ok=False)]
    lnh = an.lnh_CGM if math.isfinite(an.lnh_CGM) else an.x_min + 5.0 / 6.0
    for m in methods:
        try:
            if m == "M0":
                r = B.solve_m0(pot, an.x_min)
            elif m in ("M1", "M2"):
                # bracket search restricted to the admissible interval G < 0 (between barrier and true vacuum)
                iv = B.admissible_interval(pot, an.x_min)
                lo, hi = iv if iv else (lnh - 3, lnh + 3)
                hi = min(hi, an.x_min + 8.0)
                r = (B.solve_m1 if m == "M1" else B.solve_m2)(pot, lo, hi, n_scan=25)
            elif m == "M3":
                r = B.solve_m3(pot, an.lambda_min, lnh, max_nodes=50000)   # explicit budget; failure if exceeded
            elif m == "M4":
                # adaptive resolution: double N until both identities hold to 1e-9 (same rule as the C++ solver)
                r = None
                for N in (256, 512, 1024):
                    r = B.solve_m4(pot, an.lambda_min, lnh, N=N)
                    if r.ok and abs(r.derrick) < 1e-9 and abs(r.dilatation) < 1e-9:
                        break
        except Exception as e:  # record failures instead of hiding them
            r = B.BounceResult(m, False, info=dict(reason=repr(e)))
        row = dict(base, method=m, ok=r.ok, S=r.S, derrick=r.derrick, dilatation=r.dilatation, lnphi0=r.lnphi0,
                   lnR=r.lnR, seconds=r.seconds,
                   info=";".join(f"{k}={v}" for k, v in r.info.items() if k != "scan"))
        rows.append(row)
        append(row)
        print(f"  {name:28s} {m}: ok={r.ok!s:5} S={r.S:.12g}  derrick={r.derrick:.1e} dil={r.dilatation:.1e} "
              f"t={r.seconds:.1f}s", flush=True)
    return rows


def main():
    rows = []
    l0 = -0.0125
    # B1 constant lambda (exact S known); shooting is degenerate (every phi0 is a solution) and omitted
    rows += run_all("B1 const", P.constant_model(l0), (30, 50), methods=("M0", "M3", "M4"))
    # B2/B3 parabolic lambda: thick (q >> 1) to thin wall (q << 1); q = |l0| / l2, l2 = 5e-5 (SM-like)
    for q in ([250, 25, 2.5, 0.5] if QUICK else [250, 25, 2.5, 1.0, 0.5, 0.25, 0.1]):
        l2 = 5e-5
        rows += run_all(f"B2/3 parabola q={q}", P.parabola_model(-q * l2, l2, 39.0), (25, 55))
    # B4 SM points in the controlled region (prescription cw1), plus the tree prescription at the centre
    pts = [(125.15, 173.34), (120.0, 172.0), (130.0, 178.0), (125.15, 178.5)]
    for Mh, Mt in pts:
        rows += run_all(f"B4 SM cw1 ({Mh},{Mt})", P.sm_potential(sm.Inputs(Mh, Mt), "cw1"), None)
    rows += run_all("B4 SM tree (125.15,173.34)", P.sm_potential(sm.Inputs(125.15, 173.34), "tree"), None)
    Mtc = critical_Mt(125.15)
    print(f"critical M_t (cw1, M_h=125.15) = {Mtc:.6f} GeV", flush=True)
    for d in ([1.0, 0.1, 0.01] if QUICK else [1.0, 0.3, 0.1, 0.03, 0.01, 0.003]):
        rows += run_all(f"B4 SM cw1 Mt=Mtc+{d}", P.sm_potential(sm.Inputs(125.15, Mtc + d), "cw1"), None,
                        notes=f"Mtc={Mtc:.6f};")
    # B5 exploratory points (outside the controlled region)
    for Mh, Mt in ([(60.0, 150.0)] if QUICK else [(60.0, 150.0), (100.0, 190.0)]):
        rows += run_all(f"B5 SM cw1 ({Mh},{Mt})", P.sm_potential(sm.Inputs(Mh, Mt), "cw1"), None,
                        notes="exploratory;")
    print("wrote", os.path.join(OUT, "method_comparison.csv"))


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"total {time.time() - t0:.0f} s")
