"""Reproducibility check: recompute a seeded, stratified sample of stored results from the committed source.

1. Scan points (C++): for each committed scan, up to --per-status points per status class plus reduced-accuracy and
   Planck-sensitive points are recomputed with the given smvac_point binary and compared with the stored values
   (status, S_num, S_FL, S_CGM; relative tolerance 1e-12, the printed precision is 15 significant digits).
2. Shooting fallback (Python): --fallback points of results/controlled_cw1_fallback.csv per fallback status are
   re-solved with validation/fallback_shooting.py (status and S_fallback, relative tolerance 1e-10).
3. Near-critical model (Python): --model points of validation/results/near_critical_model.csv are recomputed
   (S_num and S_CGM2, relative tolerance 1e-9; the variational optimiser is the least reproducible element).

Usage: python validation/reproduce_sample.py --binary build_clean/smvac_point [--per-status 12] [--fallback 2]
       [--model 2]
Output: validation/results/reproduce_sample.csv and a summary on stdout. Exit status 1 if any mismatch.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCANS = {"controlled_cw1": ["--prescription", "cw1", "--matching", "buttazzo2013"],
         "controlled_tree": ["--prescription", "tree", "--matching", "buttazzo2013"],
         "exploratory_cw1": ["--prescription", "cw1", "--matching", "extended"]}


def load(name):
    p = os.path.join(ROOT, "results", name + ".csv")
    return pd.read_csv(p if os.path.exists(p) else p + ".gz", comment="#")


def close(a, b, rel=1e-12):
    """a: recomputed value (None = JSON null, which encodes any non-finite number); b: stored value."""
    b = float(b)
    a_finite = a is not None and np.isfinite(a)
    if not a_finite or not np.isfinite(b):
        return not a_finite and not np.isfinite(b)
    return abs(a - b) <= rel * max(abs(b), 1e-300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--per-status", type=int, default=12)
    ap.add_argument("--fallback", type=int, default=2)
    ap.add_argument("--model", type=int, default=2)
    a = ap.parse_args()
    binary = os.path.abspath(a.binary + (".exe" if os.name == "nt" and not a.binary.endswith(".exe") else ""))
    rng = np.random.default_rng(20261109)
    rows, bad = [], 0
    for name, opts in SCANS.items():
        df = load(name)
        groups = [g for _, g in df.groupby("status")] + [df[df.reduced_accuracy == 1], df[df.planck_sensitive == 1]]
        sample = pd.concat([g.iloc[rng.choice(len(g), size=min(a.per_status, len(g)), replace=False)]
                            for g in groups if len(g)]).drop_duplicates(["Mh", "Mt"])
        for r in sample.itertuples():
            j = json.loads(subprocess.run([binary, repr(float(r.Mh)), repr(float(r.Mt)), *opts], capture_output=True,
                                          text=True, check=True).stdout)
            ok = (j["status"] == r.status and close(j["S_numerical"], r.S_num) and close(j["S_FL_lower_bound"], r.S_FL)
                  and close(j["S_CGM_upper_estimate"], r.S_CGM))
            bad += 0 if ok else 1
            rows.append(dict(scan=name, Mh=r.Mh, Mt=r.Mt, status_committed=r.status, status_recomputed=j["status"],
                             S_committed=r.S_num, S_recomputed=j["S_numerical"], match=ok))
    # 2. shooting fallback
    sys.path.insert(0, os.path.dirname(__file__))
    fbp = os.path.join(ROOT, "results", "controlled_cw1_fallback.csv")
    if a.fallback > 0 and os.path.exists(fbp):
        from fallback_shooting import solve_point
        fb = pd.read_csv(fbp, comment="#")
        for st, g in fb.groupby("fallback_status"):
            for r in g.iloc[rng.choice(len(g), size=min(a.fallback, len(g)), replace=False)].itertuples():
                new = solve_point((float(r.Mh), float(r.Mt), "cw1", "buttazzo2013"))
                ok = new["fallback_status"] == st and (st not in ("metastable", "unstable")
                                                       or close(new["S_fallback"], r.S_fallback, 1e-10))
                bad += 0 if ok else 1
                rows.append(dict(scan="fallback_cw1", Mh=r.Mh, Mt=r.Mt, status_committed=st,
                                 status_recomputed=new["fallback_status"], S_committed=getattr(r, "S_fallback", np.nan),
                                 S_recomputed=new.get("S_fallback", np.nan), match=ok))
    # 3. near-critical model family
    mp = os.path.join(ROOT, "validation", "results", "near_critical_model.csv")
    if a.model > 0 and os.path.exists(mp):
        from near_critical import model_point
        m = pd.read_csv(mp)
        for r in m.iloc[rng.choice(len(m), size=min(a.model, len(m)), replace=False)].itertuples():
            new = model_point(float(r.q))
            ok = close(new["S_num"], r.S_num, 1e-9) and (
                (not np.isfinite(r.S_CGM2) and not np.isfinite(new["S_CGM2"])) or close(new["S_CGM2"], r.S_CGM2, 1e-9))
            bad += 0 if ok else 1
            rows.append(dict(scan="near_critical_model", Mh=np.nan, Mt=np.nan, status_committed=f"q={r.q:g}",
                             status_recomputed=f"q={r.q:g}", S_committed=r.S_num, S_recomputed=new["S_num"], match=ok))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(ROOT, "validation", "results", "reproduce_sample.csv"), index=False)
    print(out.groupby("scan").match.agg(["count", "sum"]).rename(columns={"sum": "matched"}))
    print("mismatches:", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
