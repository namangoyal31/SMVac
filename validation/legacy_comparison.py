"""Evaluate the legacy solver on a coarse version of the controlled grid and compare with the new results.

Requires the legacy point driver built as described in audit/README.md (path given with --driver).
Output: validation/results/legacy_comparison.csv with legacy status/S_exact/S_approx next to the new
status/S_num/S_FL (prescription cw1) at the same (M_h, M_t).
Usage: python validation/legacy_comparison.py --driver path/to/legacy_driver[.exe]
"""
from __future__ import annotations

import argparse
import io
import os
import subprocess

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LEGACY_STATUS = {1: "stable", 2: "metastable", 3: "unstable", 4: "non-perturbative"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", required=True)
    ap.add_argument("--step", type=float, default=0.5)
    a = ap.parse_args()
    path = os.path.join(ROOT, "results", "controlled_cw1.csv")
    new = pd.read_csv(path if os.path.exists(path) else path + ".gz", comment="#")
    mh = np.round(np.arange(120.0, 130.0 + 1e-9, a.step), 6)
    mt = np.round(np.arange(166.0, 182.0 + 1e-9, a.step), 6)
    pts = [(h, t) for t in mt for h in mh]
    inp = "".join(f"{float(h)!r} {float(t)!r}\n" for h, t in pts)   # plain floats (numpy 2 repr differs)
    out = subprocess.run([os.path.abspath(a.driver)], input=inp, capture_output=True, text=True, check=True).stdout
    leg = pd.read_csv(io.StringIO(out))
    leg["legacy_status"] = leg.Stability.map(LEGACY_STATUS)
    leg = leg.rename(columns={"S_exact": "legacy_S_exact", "S_approx": "legacy_S_approx"})
    new["Mh_r"], new["Mt_r"] = new.Mh.round(6), new.Mt.round(6)
    leg["Mh_r"], leg["Mt_r"] = leg.Mh.round(6), leg.Mt.round(6)
    m = leg.merge(new, on=["Mh_r", "Mt_r"], suffixes=("_leg", ""))
    cols = ["Mh", "Mt", "legacy_status", "legacy_S_exact", "legacy_S_approx", "status", "S_num", "S_FL", "lambda_min"]
    m[cols].to_csv(os.path.join(ROOT, "validation", "results", "legacy_comparison.csv"), index=False)
    print(pd.crosstab(m.legacy_status, m.status))


if __name__ == "__main__":
    main()
