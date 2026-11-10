"""Extract the quantitative results quoted in README.md and docs/VALIDATION_REPORT.md from the scan outputs.

Writes results/summary.json and prints a human-readable digest. Every number quoted in the documentation should be
traceable to this script.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# H_0 = 67.4 km/s/Mpc in GeV (Buttazzo et al. eq. 82 convention, as in src/decay.cpp). With Lambda_B <= M_Pl,
# p = 0.15 Lambda_B^4 e^-S / H_0^4 < 1 whenever S > S_CRIT_PLANCK; since S >= S_FL for every finite-action bounce
# (rigorous bound), S_FL > S_CRIT_PLANCK classifies a point as metastable without solving the bounce equation.
H0_GEV = 67.4 / 3.0856775814913673e19 * 6.582119569e-25
S_CRIT_PLANCK = float(np.log(0.15) + 4.0 * np.log(1.22e19 / H0_GEV))


def load(name):
    p = os.path.join(ROOT, "results", name)
    if not os.path.exists(p) and os.path.exists(p + ".gz"):
        p += ".gz"                     # committed scan outputs are gzip-compressed
    return pd.read_csv(p, comment="#") if os.path.exists(p) else None


def quantiles(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return None
    return {"n": int(len(x)), "min": float(x.min()), "median": float(np.median(x)), "p90": float(np.quantile(x, 0.9)),
            "max": float(x.max())}


def boundary_along(df, Mh, col, level=0.0):
    """M_t at which column `col` crosses `level` along fixed M_h (linear interpolation between grid nodes)."""
    d = df[np.isclose(df.Mh, Mh)].sort_values("Mt")
    v, t = d[col].values, d.Mt.values
    out = []
    for i in range(len(v) - 1):
        if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and (v[i] - level) * (v[i + 1] - level) < 0:
            out.append(float(t[i] + (level - v[i]) * (t[i + 1] - t[i]) / (v[i + 1] - v[i])))
    return out


def analyze(df, label):
    s = {"label": label, "points": int(len(df)), "status_counts": df.status.value_counts().to_dict()}
    b = df[df.S_num.notna()].copy()
    clean = b[(b.planck_sensitive == 0) & (b.reduced_accuracy == 0)]
    s["bounces"] = int(len(b))
    s["planck_sensitive"] = int((df.planck_sensitive == 1).sum())
    s["reduced_accuracy"] = int((b.reduced_accuracy == 1).sum())
    s["max_rel_uncertainty_clean"] = float((clean.dS_num / clean.S_num).max()) if len(clean) else None
    s["max_identity_residual_clean"] = float(np.maximum(clean.derrick.abs(), clean.dilatation.abs()).max()) if len(clean) else None
    s["bound_violations"] = int(((b.S_FL > b.S_num * (1 + 1e-12)) | (b.S_num > b.S_CGM * (1 + 1e-9))).sum())
    for col in ("S_FL", "S_curv", "S_CGM"):
        s[f"rel_error_{col}"] = quantiles(((b[col] - b.S_num) / b.S_num).abs())
        s[f"rel_error_{col}_clean"] = quantiles(((clean[col] - clean.S_num) / clean.S_num).abs())
    # error of S_FL in the deep-metastable part (S_num < 5000) vs near-critical part
    s["rel_error_S_FL_S_lt_5000"] = quantiles(((b.S_num - b.S_FL) / b.S_num)[b.S_num < 5000])
    s["rel_error_S_FL_S_gt_1e4"] = quantiles(((b.S_num - b.S_FL) / b.S_num)[b.S_num > 1e4])
    # absolute action shift relevant for the decay exponent near the metastability boundary (S ~ 500-600)
    nb = b[(b.S_num > 450) & (b.S_num < 650)]
    s["abs_error_S_FL_near_metastability_boundary"] = quantiles(nb.S_num - nb.S_FL)
    s["threshold_scale_shift"] = quantiles((b.log10_p_FL_lit - b.log10_p_FL_bscale) * np.log(10))
    for Mh in (125.15, 125.0, 125.2):
        if np.isclose(df.Mh, Mh).any():
            s[f"boundaries_at_Mh_{Mh}"] = {
                "stability_lambda_min_0": boundary_along(df, Mh, "lambda_min"),
                "metastability_numerical": boundary_along(df, Mh, "log10_p_num"),
                "metastability_FL_literature_scale": boundary_along(df, Mh, "log10_p_FL_lit"),
                "metastability_FL_bounce_scale": boundary_along(df, Mh, "log10_p_FL_bscale"),
            }
    return s


def main():
    out = {}
    cw1, tree, expl = load("controlled_cw1.csv"), load("controlled_tree.csv"), load("exploratory_cw1.csv")
    if cw1 is not None:
        out["controlled_cw1"] = analyze(cw1, "controlled region, cw1")
    if tree is not None:
        out["controlled_tree"] = analyze(tree, "controlled region, tree")
    if cw1 is not None and tree is not None:
        m = cw1.merge(tree, on=["Mh", "Mt"], suffixes=("_cw1", "_tree"))
        mm = m[m.S_num_cw1.notna() & m.S_num_tree.notna()]
        out["prescription_systematic"] = {
            "rel_diff_S_num_tree_vs_cw1": quantiles(((mm.S_num_tree - mm.S_num_cw1) / mm.S_num_cw1).abs()),
            "status_disagreements": int((m.status_cw1 != m.status_tree).sum()),
            "central_point": {k: float(m[(np.isclose(m.Mh, 125.2)) & (np.isclose(m.Mt, 173.3))][k].iloc[0])
                              for k in ("S_num_cw1", "S_num_tree")} if len(m[(np.isclose(m.Mh, 125.2)) & (np.isclose(m.Mt, 173.3))]) else None,
        }
    if expl is not None:
        out["exploratory_cw1"] = {"points": int(len(expl)), "status_counts": expl.status.value_counts().to_dict(),
                                  "planck_sensitive": int(expl.planck_sensitive.sum()),
                                  "reduced_accuracy": int(expl.reduced_accuracy.sum())}
    fb_path = os.path.join(ROOT, "results", "controlled_cw1_fallback.csv")
    if os.path.exists(fb_path):
        fb = pd.read_csv(fb_path, comment="#")
        acc = fb[fb.fallback_status.isin(["metastable", "unstable"])]
        out["fallback_cw1"] = {
            "points": int(len(fb)), "fallback_status_counts": fb.fallback_status.value_counts().to_dict(),
            "max_rel_diff_M2_M1": float(((acc.S_M2 - acc.S_M1).abs() / acc.S_M2).max()) if len(acc) else None,
            "max_abs_derrick": float(np.maximum(acc.derrick_M2.abs(), acc.derrick_M1.abs()).max()) if len(acc) else None,
            "planck_sensitive": int(acc.planck_sensitive.sum()) if len(acc) else 0,
            "log10_phi0_range": [float(acc.log10_phi0.min()), float(acc.log10_phi0.max())] if len(acc) else None,
            "rel_error_S_FL": quantiles((acc.S_fallback - acc.S_FL) / acc.S_fallback),
        }
        nfb = fb[fb.fallback_status == "no-finite-bounce"]
        out["fallback_cw1"]["no_finite_bounce"] = {
            "points": int(len(nfb)),
            "Mh_range": [float(nfb.Mh.min()), float(nfb.Mh.max())], "Mt_range": [float(nfb.Mt.min()), float(nfb.Mt.max())],
            "log10_mu_min_range": [float(nfb.log10_mu_min.min()), float(nfb.log10_mu_min.max())],
            "plateau_rise_max": float(nfb.plateau_rise.max()),
            "S_FL_min": float(nfb.S_FL.min()),
            "S_plateau_infimum_over_S_FL_max": float((nfb.S_plateau_infimum / nfb.S_FL).max()),
        }
        if cw1 is not None:
            # final classification of the controlled cw1 scan after the shooting fallback, and the points without a
            # finite bounce (lambda minimum on, or within ~0.2 decades below, the frozen boundary at M_Pl)
            merged = cw1.status.copy()
            key = {(round(a, 6), round(b, 6)): st for a, b, st in zip(fb.Mh, fb.Mt, fb.fallback_status)}
            for j in cw1.index[cw1.status == "solver-failure"]:
                merged[j] = key.get((round(cw1.Mh[j], 6), round(cw1.Mt[j], 6)), "solver-failure")
            out["controlled_cw1"]["final_status_counts_with_fallback"] = merged.value_counts().to_dict()
            nob = cw1[cw1.status == "no-bounce"]
            plateau = pd.concat([nob[["Mh", "Mt", "S_FL", "log10_mu_min"]], nfb[["Mh", "Mt", "S_FL", "log10_mu_min"]]])
            out["controlled_cw1"]["planck_plateau_points"] = {
                "no_bounce_minimum_at_M_Pl": int(len(nob)), "no_finite_bounce_minimum_below_M_Pl": int(len(nfb)),
                "Mh_min": float(plateau.Mh.min()), "S_FL_min": float(plateau.S_FL.min()),
                "S_threshold_for_p_lt_1_with_LambdaB_le_M_Pl": S_CRIT_PLANCK,
                "all_metastable_by_rigorous_bound": bool(plateau.S_FL.min() > S_CRIT_PLANCK),
            }
    with open(os.path.join(ROOT, "results", "summary.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
