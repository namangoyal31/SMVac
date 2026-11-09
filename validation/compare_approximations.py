"""Consistent comparison of analytical approximations with the numerical reference (post-processing only).

Reads the committed datasets (no solver is run):
  results/controlled_cw1.csv.gz, results/controlled_tree.csv.gz, results/exploratory_cw1.csv.gz,
  validation/results/near_critical_sm.csv, validation/results/near_critical_model.csv
and writes validation/results/approximation_comparison.json and approximation_comparison.md.

Approximations, all evaluated from stored columns (S_FL = 8 pi^2/(3|lambda_min|), q = |lambda_min|/lambda''_min):
  FL    S_FL                                   leading Fubini-Lipatov formula; rigorous lower bound
  E1    S_FL (1 + c1/q),  c1 = 13/72           first-order formula of J.R. Espinosa, JCAP 06 (2020) 052, eq. (3.21),
                                               with lambda_0 = |lambda_min|, beta'_lambda = lambda''_min
  E2    S_FL (1 + c1/q + c2/q^2)               second-order term of Espinosa eq. (6.5) (kappa = 0, L_0 = 5/6); exact only
                                               for lambda quadratic in ln(phi) (supplementary column)
  curv  S_FL / (1 - c1/q)                      resummed first-order form ("curvature correction", stored column S_curv)
  CGM   8 pi^2 / (3 |min_h lambdabar(h)|)      one-parameter (Fubini-shape) variational estimate (stored column)
  CGM2  two-parameter variational estimate     near-critical datasets only (stored column)
Error metric: signed relative error e = (S_est - S_num)/S_num and its absolute value |e|. A non-finite or missing
estimate at a solved point (the variational/resummed estimates diverge) is counted as |e| = infinity.

Selections (identical for every approximation within a dataset):
  scans: S_num finite and > 0, lambda'' > 0, not Planck-sensitive, not reduced-accuracy ("standard tier");
  near-critical SM lines: solved points (S_num finite), per line; model family: all 17 points.
Regions: q >= 10 (bulk / thick wall), 1 <= q < 10, 0.1 <= q < 1, q < 0.1 (near-critical), and "all".
The exploratory scan is reported as its own region (uncontrolled extrapolation of the matching).
Usage: python validation/compare_approximations.py
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(ROOT, "validation", "results")
C1 = 13.0 / 72.0
BANDS = [("q>=10", 10.0, math.inf), ("1<=q<10", 1.0, 10.0), ("0.1<=q<1", 0.1, 1.0), ("q<0.1", 0.0, 0.1),
         ("all", 0.0, math.inf)]


def espinosa_second_order_coefficient() -> float:
    """c2 such that S/S_FL = 1 + c1/q + c2/q^2, from Espinosa (2020) eq. (6.5) at kappa = 0 and L_0 = 5/6.

    In eq. (6.5) the O(beta'^2) term is (pi^2 beta'^2 / (9 lambda_0^3)) B(L_0); dividing by S_FL = 8 pi^2/(3 lambda_0)
    gives B/24 with q = lambda_0/beta'. The O(beta') shift of L_0 contributes only at O(beta'^3) (S is stationary
    in L_0), so L_0 = 5/6 may be used."""
    from scipy.integrate import quad
    from scipy.special import spence
    li2 = lambda x: spence(1.0 - x)          # dilogarithm Li_2(x)
    k0 = quad(lambda x: (1 - 6 * x + 6 * x * x) * li2(x) ** 2, 0.0, 1.0, epsabs=0.0, epsrel=1e-13)[0]
    L = 5.0 / 6.0
    b = (18 * k0 + 553 / 36 - 11 * math.pi ** 2 / 6 + 2 * (math.pi ** 2 - 12) * L + 8 * L ** 2 - 8 * L ** 3
         + 6 * L ** 4)
    return b / 24.0


def load(path):
    p = os.path.join(ROOT, path)
    if not os.path.exists(p) and os.path.exists(p + ".gz"):
        p += ".gz"
    return pd.read_csv(p, comment="#")


def estimates(df, c2):
    """Dictionary name -> array of estimates (non-finite where the estimate does not exist)."""
    q = (df.lambda_min.abs() / df.lambda2).to_numpy(float)
    sfl = df.S_FL.to_numpy(float)
    clean = lambda a: np.where(np.isfinite(a) & (a > 0), a, np.inf)
    out = {"FL": sfl, "E1": sfl * (1 + C1 / q), "E2": sfl * (1 + C1 / q + c2 / q ** 2),
           "curv": clean(df.S_curv.to_numpy(float)), "CGM": clean(df.S_CGM.to_numpy(float))}
    if "S_CGM2" in df:
        out["CGM2"] = clean(df.S_CGM2.to_numpy(float))
    return q, out


def stats(est, snum, q, meta):
    """Error statistics for one approximation on one selection (arrays of equal length)."""
    n = int(len(snum))
    if n == 0:
        return {"n": 0}
    e = (est - snum) / snum
    a = np.abs(e)
    fin = np.isfinite(a)
    d = {"n": n, "n_nonfinite": int((~fin).sum()),
         # empirical quantiles without interpolation, so that infinite errors do not produce inf - inf = nan
         "median_abs": float(np.quantile(a, 0.5, method="inverted_cdf")),
         "p90_abs": float(np.quantile(a, 0.9, method="inverted_cdf")), "max_abs": float(a.max()),
         "median_signed": float(np.quantile(e, 0.5, method="inverted_cdf"))}
    if fin.any():
        i = int(np.argmax(np.where(fin, a, -1.0)))
        d.update(max_abs_finite=float(a[i]), signed_at_max=float(e[i]), q_at_max=float(q[i]),
                 S_num_at_max=float(snum[i]), signed_min=float(e[fin].min()), signed_max=float(e[fin].max()))
        d.update({k: float(v[i]) for k, v in meta.items()})
    if (~fin).any():
        d["q_max_nonfinite"] = float(q[~fin].max())
    return d


def compare(df, c2, meta_cols):
    q, est = estimates(df, c2)
    snum = df.S_num.to_numpy(float)
    out = {}
    for band, lo, hi in BANDS:
        m = (q >= lo) & (q < hi)
        meta = {c: df[c].to_numpy(float)[m] for c in meta_cols}
        out[band] = {name: stats(v[m], snum[m], q[m], meta) for name, v in est.items()}
    return out


def standard_tier(df):
    return df[np.isfinite(df.S_num) & (df.S_num > 0) & (df.lambda2 > 0) & (df.planck_sensitive == 0)
              & (df.reduced_accuracy == 0)].reset_index(drop=True)


def coverage(raw, strict, fallback=None):
    """Point accounting for one dataset: why each row is or is not in the comparison.
    strict=True: standard tier (scans); strict=False: all solved points with lambda'' > 0 (near-critical lines)."""
    solved = np.isfinite(raw.S_num) & (raw.S_num > 0)
    lam_ok = solved & (raw.lambda2 > 0)
    planck = lam_ok & (raw.planck_sensitive.astype(int) == 1)
    reduced = lam_ok & ~planck & (raw.reduced_accuracy.astype(int) == 1)
    included = (lam_ok & ~planck & ~reduced) if strict else lam_ok
    c = {"rows": int(len(raw)), "status_counts": {str(k): int(v) for k, v in raw.status.value_counts().items()},
         "no_numerical_action": int((~solved).sum()), "solved": int(solved.sum()),
         "solved_but_lambda2_not_positive": int((solved & ~lam_ok).sum()),
         "solved_planck_sensitive": int(planck.sum()), "solved_reduced_accuracy_not_planck": int(reduced.sum()),
         "included": int(included.sum()),
         "inclusion_rule": "standard tier" if strict else "all solved points with lambda'' > 0 (Planck-sensitive "
                                                          "and reduced-accuracy points included)"}
    if fallback is not None:          # classification of collocation failures by the shooting fallback (not compared)
        c["solver_failures_fallback_status"] = {str(k): int(v) for k, v in fallback.fallback_status.value_counts().items()}
    return c


def run():
    c2 = espinosa_second_order_coefficient()
    res = {"definitions": {"c1": C1, "c2": c2, "error": "e = (S_est - S_num)/S_num; non-finite estimate -> |e| = inf",
                           "standard_tier": "S_num finite > 0, lambda'' > 0, planck_sensitive = 0, reduced_accuracy = 0"}}
    res["coverage"] = {}
    for name in ("controlled_cw1", "controlled_tree", "exploratory_cw1"):
        raw = load(f"results/{name}.csv")
        fb = load("results/controlled_cw1_fallback.csv") if name == "controlled_cw1" else None
        res["coverage"][name] = coverage(raw, True, fb)
        d = standard_tier(raw)
        res[name] = compare(d, c2, ["Mh", "Mt"])
    both = pd.concat([standard_tier(load("results/controlled_cw1.csv")), standard_tier(load("results/controlled_tree.csv"))],
                     ignore_index=True)
    res["controlled_both"] = compare(both, c2, ["Mh", "Mt"])
    sm = load("validation/results/near_critical_sm.csv")
    res["coverage"]["near_critical_sm"] = coverage(sm, False)
    sm = sm[np.isfinite(sm.S_num) & (sm.S_num > 0) & (sm.lambda2 > 0)]
    res["near_critical_sm"] = {f"{p} Mh={mh:g}": compare(g.reset_index(drop=True), c2, ["Mh", "Mt", "delta"])
                               for (p, mh), g in sm.groupby(["prescription", "Mh"])}
    res["near_critical_sm_all"] = compare(sm.reset_index(drop=True), c2, ["Mh", "Mt", "delta"])
    model = load("validation/results/near_critical_model.csv")
    model = model.assign(lambda_min=-model.q * model.l2, lambda2=model.l2)
    res["coverage"]["model_family"] = {"rows": int(len(model)), "included": int(np.isfinite(model.S_num).sum()),
                                       "inclusion_rule": "all model points (idealised potential; no Planck scale)"}
    res["model_family"] = compare(model.reset_index(drop=True), c2, [])
    q, est = estimates(model, c2)
    res["model_family_table"] = [{"q": float(qq), **{k: float((v[i] - model.S_num[i]) / model.S_num[i])
                                                    for k, v in est.items()}} for i, qq in enumerate(q)]
    return res


def fmt(x):
    if x is None:
        return "–"
    if not math.isfinite(x):
        return "∞"
    return f"{x:.1e}" if abs(x) < 1e-2 else f"{x:.3f}"


def markdown(res):
    names = ["FL", "E1", "E2", "curv", "CGM", "CGM2"]
    lines = ["# Approximation comparison (generated by validation/compare_approximations.py)", "",
             f"c1 = 13/72 = {res['definitions']['c1']:.6f}; c2 (Espinosa 2020 eq. 6.5) = {res['definitions']['c2']:.6f}.",
             "Entries: median |e| / max |e| (∞: the estimate is infinite at some point; the largest finite error and the",
             "count of infinite cases are in the JSON file). e = (S_est - S_num)/S_num.", "",
             "## Coverage (which points enter each comparison)", "",
             "| dataset | rows | no numerical action | solved | Planck-sensitive (solved) | reduced accuracy (solved) | included | rule |",
             "|---|---|---|---|---|---|---|---|"]
    for ds, c in res["coverage"].items():
        lines.append(f"| {ds} | {c['rows']} | {c.get('no_numerical_action', 0)} | {c.get('solved', c['rows'])} | "
                     f"{c.get('solved_planck_sensitive', 0)} | {c.get('solved_reduced_accuracy_not_planck', 0)} | "
                     f"{c['included']} | {c['inclusion_rule']} |")
    lines += ["", "Rows without a numerical action are stable points, points without a finite bounce, collocation "
              "failures, non-perturbative and invalid inputs (status counts in the JSON file). The 360 `cw1` collocation "
              "failures are classified by the shooting fallback (78 metastable, all Planck-sensitive; 282 no finite "
              "bounce); none of them enters the comparison.", ""]
    for ds in ("controlled_both", "controlled_cw1", "controlled_tree", "exploratory_cw1", "near_critical_sm_all",
               "model_family"):
        lines += [f"## {ds}", "", "| region | n | " + " | ".join(names) + " |", "|---" * (len(names) + 2) + "|"]
        for band, _, _ in BANDS:
            row = res[ds][band]
            n = row["FL"]["n"]
            if n == 0:
                continue
            cells = []
            for k in names:
                s = row.get(k)
                cells.append("–" if s is None else f"{fmt(s['median_abs'])} / {fmt(s['max_abs'])}")
            lines.append(f"| {band} | {n} | " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines) + "\n"


def main():
    res = run()
    with open(os.path.join(RES, "approximation_comparison.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, default=float)
    md = markdown(res)
    with open(os.path.join(RES, "approximation_comparison.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(md)
    print(md)


if __name__ == "__main__":
    main()
