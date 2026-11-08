"""Quantitative study of the approach to the absolute-stability boundary (thin-wall regime).

Part A, model family lambda(x) = l0 + (l2/2)(x - x0)^2 with q = |l0|/l2 from 1000 down to 0.05:
  numerical reference (M4 with adaptive N, cross-checked by M2 shooting), S_FL, S_curv, S_CGM and the two-parameter
  variational estimate S_CGM2 (python/smvac/reference/analytic.py: cgm2_estimate).
  For this family S_curv = S_CGM = 8 pi^2 / (3 l2 (q - 13/72)) exactly, which diverges at q = 13/72; S_CGM2 is
  finite only for q > 3/32 (python/smvac/reference/analytic.py, two-parameter section).
Part B, SM lines approaching M_t^crit (min lambda_P = 0) at fixed M_h:
  cw1 at M_h = 120, 125.15, 128 GeV and tree at M_h = 125.15 GeV, delta = M_t - M_t^crit from 1e-4 to 3 GeV.
  Numerical reference from the C++ binary (build/smvac_point) with the shooting fallback where it fails; estimates
  from the C++ output; S_CGM2 from the Python reference.
Outputs: validation/results/near_critical_model.csv, validation/results/near_critical_sm.csv
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.optimize import brentq

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, os.path.dirname(__file__))
BIN = os.path.join(ROOT, "build", "smvac_point" + (".exe" if os.name == "nt" else ""))
RES = os.path.join(ROOT, "validation", "results")


def model_point(q):
    from smvac.reference import analytic as A, bounce as B, potential as P
    l2 = 5e-5
    pot = P.parabola_model(-q * l2, l2, 39.0)
    an = A.estimates(pot, (15, 65))
    lnh = an.lnh_CGM if math.isfinite(an.lnh_CGM) else an.x_min + 5.0 / 6.0
    r4 = None
    for N in (256, 512, 1024, 2048):
        r4 = B.solve_m4(pot, an.lambda_min, lnh, N=N)
        if r4.ok and abs(r4.derrick) < 1e-9 and abs(r4.dilatation) < 1e-9:
            break
    iv = B.admissible_interval(pot, an.x_min, half=10.0)
    r2 = B.solve_m2(pot, iv[0], iv[1], n_scan=41)
    S2, lnh2, n2 = A.cgm2_estimate(pot, lnh)
    ok4 = r4.ok and abs(r4.derrick) < 1e-9 and abs(r4.dilatation) < 1e-9
    S = r4.S if ok4 else (r2.S if r2.ok else float("nan"))
    return dict(q=q, l2=l2, S_num=S, S_M4=r4.S if ok4 else float("nan"), N_M4=r4.info.get("N"), S_M2=r2.S if r2.ok else float("nan"),
                S_FL=an.S_FL, S_curv=an.S_curv, S_CGM=an.S_CGM, S_CGM2=S2, n_CGM2=n2,
                exact_CGM=8 * math.pi ** 2 / (3 * l2 * (q - 13 / 72)) if q > 13 / 72 else float("inf"))


def mt_crit(Mh, presc):
    from smvac.reference import analytic as A, potential as P, sm
    return brentq(lambda Mt: A.global_min(P.sm_potential(sm.Inputs(Mh, Mt), presc))[0], 166.0, 178.0, xtol=1e-8)


def sm_point(args):
    Mh, Mt, Mtc, presc = args
    from smvac.reference import analytic as A, potential as P, sm
    j = json.loads(subprocess.run([BIN, repr(Mh), repr(Mt), "--prescription", presc], capture_output=True, text=True,
                                  check=True).stdout)
    row = dict(Mh=Mh, Mt=Mt, delta=Mt - Mtc, prescription=presc, status=j["status"], method="chebyshev-collocation",
               S_num=j["S_numerical"], dS=j["dS_numerical"], S_FL=j["S_FL_lower_bound"], S_curv=j["S_curv"],
               S_CGM=j["S_CGM_upper_estimate"], lambda_min=j["lambda_min"], lambda2=j["lambda2_at_min"],
               log10_mu_min=j["log10_mu_lambda_min"], log10_phi0=j["log10_phi0"], planck_sensitive=j["planck_sensitive"],
               reduced_accuracy=j["reduced_accuracy"], message=j["message"])
    if row["status"] == "solver-failure":
        from fallback_shooting import solve_point
        fb = solve_point((Mh, Mt, presc, "buttazzo2013"))
        if fb.get("fallback_status") in ("metastable", "unstable"):
            row.update(status=fb["fallback_status"], method="shooting-fallback", S_num=fb["S_fallback"],
                       dS=fb["dS_fallback"], log10_phi0=fb["log10_phi0"], planck_sensitive=bool(fb["planck_sensitive"]))
        else:
            row.update(status=fb.get("fallback_status", "failed"), method="shooting-fallback (no solution)")
    pot = P.sm_potential(sm.Inputs(Mh, Mt), presc)
    an = A.estimates(pot)
    lnh = an.lnh_CGM if math.isfinite(an.lnh_CGM) else an.x_min + 5.0 / 6.0
    row["S_CGM2"], _, row["n_CGM2"] = A.cgm2_estimate(pot, lnh)
    return row


def main():
    workers = max(1, int(os.environ.get("SMVAC_WORKERS", (os.cpu_count() or 2) - 2)))
    qs = [1000.0, 300.0, 100.0, 30.0, 10.0, 3.0, 1.0, 0.5, 0.35, 0.25, 0.2, 0.15, 0.12, 0.1, 0.085, 0.07, 0.05]
    jobs = []
    for Mh, presc in ((120.0, "cw1"), (125.15, "cw1"), (128.0, "cw1"), (125.15, "tree")):
        Mtc = mt_crit(Mh, presc)
        print(f"{presc} M_h = {Mh}: M_t^crit = {Mtc:.8f} GeV", flush=True)
        jobs += [(Mh, Mtc + float(d), Mtc, presc) for d in np.logspace(-4, math.log10(3.0), 17)]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        model = list(ex.map(model_point, qs, chunksize=1))
        sm_rows = list(ex.map(sm_point, jobs, chunksize=1))
    pd.DataFrame(model).to_csv(os.path.join(RES, "near_critical_model.csv"), index=False)
    pd.DataFrame(sm_rows).to_csv(os.path.join(RES, "near_critical_sm.csv"), index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(model)[["q", "S_num", "S_M2", "S_FL", "S_CGM", "S_CGM2", "n_CGM2"]].to_string(index=False))
    print(pd.DataFrame(sm_rows)[["prescription", "Mh", "delta", "status", "method", "S_num", "S_FL", "S_CGM", "S_CGM2",
                                 "log10_phi0", "planck_sensitive"]].to_string(index=False))


if __name__ == "__main__":
    main()
