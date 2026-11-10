"""Quantitative characterisation of the near-critical (thin-wall) regime.

Reads validation/results/near_critical_model.csv, near_critical_sm.csv and the controlled scans; writes
validation/results/near_critical_summary.json. Relative errors are |S_est - S_num| / S_num; for the lower bound S_FL
this is capped at 1, for the upper estimates it is unbounded, so the symmetric measure |ln(S_est/S_num)| is also given.

Model family lambda = l2 [(x - x0)^2/2 - q]:
  * thick-wall asymptotics: q (S - S_FL)/S_FL -> 13/72 (derived); q^2 x error of S_CGM and S_CGM2 (fitted, q >= 30);
    (n_opt - 1) q for the two-parameter trial;
  * thresholds q at which each estimate's error reaches 1%, 10%, 100%; crossovers q* where an estimate becomes worse
    than S_FL; the finiteness thresholds of S_CGM (13/72) and S_CGM2 (3/32), compared with the data;
  * thin-wall side: local exponent d ln S / d ln q, and the ratio S_num / S_TW to the leading-order Coleman thin-wall
    action S_TW = 27 pi^2 sigma^4 / (2 eps^3) with sigma = sqrt(l2) phi_-^3 / 18, eps = q l2 phi_-^4 / 4, i.e.
    S_TW = (27 * 64 / (2 * 18^4)) pi^2 / (l2 q^3) (degenerate potential l2 (x - x0)^2 phi^4 / 8; leading order in q).
SM lines (fixed M_h, M_t -> M_t^crit, delta = M_t - M_t^crit):
  * status counts (including Planck plateau and solver failures), the map q(delta), thresholds in delta and q, local
    scaling d ln S / d ln delta, Planck sensitivity;
  * universality: deviation of the SM errors from the model curve at the same q.
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(ROOT, "validation", "results")
ESTIMATES = ("S_FL", "S_curv", "S_CGM", "S_CGM2")
TW_COEFF = 27.0 * 64.0 / (2.0 * 18.0 ** 4) * math.pi ** 2      # S_TW = TW_COEFF / (l2 q^3)


def crossing(x, y, level):
    """First x (scanning from large to small x) at which y reaches `level`, log-log interpolation. NaN if none.
    An infinite y beyond a finite one counts as a crossing at the last finite point's neighbour (returned as the
    geometric mean of the two x values, flagged by the caller through the bracketing table)."""
    order = np.argsort(x)[::-1]
    x, y = np.asarray(x, float)[order], np.asarray(y, float)[order]
    for i in range(len(x) - 1):
        a, b = y[i], y[i + 1]
        if not (np.isfinite(a) and a > 0):
            continue
        if not np.isfinite(b):
            if a < level:
                return float(math.sqrt(x[i] * x[i + 1]))
            continue
        if b > 0 and (a - level) * (b - level) <= 0 and a != b:
            t = (math.log(level) - math.log(a)) / (math.log(b) - math.log(a))
            return float(math.exp(math.log(x[i]) + t * (math.log(x[i + 1]) - math.log(x[i]))))
    return float("nan")


def errors(df):
    """Relative and log errors; a missing (NaN) or infinite estimate at a solved point means 'no V < 0' -> inf."""
    rel, log = {}, {}
    for col in ESTIMATES:
        if col in df:
            est = df[col].where(np.isfinite(df[col]) & (df[col] > 0), np.inf)
            rel[col] = ((est - df.S_num).abs() / df.S_num).where(np.isfinite(est), np.inf)
            log[col] = np.abs(np.log(est / df.S_num)).where(np.isfinite(est), np.inf)
    return rel, log


def fnum(v):
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)


def model_summary(m):
    s = {}
    rel, log = errors(m)
    thick = m[m.q >= 30]
    s["thick_wall_q_times_FL_error_over_13_72"] = {f"{q:g}": float((r.S_num / r.S_FL - 1) * q / (13 / 72))
                                                   for q, r in zip(m.q, m.itertuples()) if q >= 10}
    s["thick_wall_coeff_q2_err_CGM"] = float(np.median(rel["S_CGM"][m.q >= 30] * m.q[m.q >= 30] ** 2))
    s["thick_wall_coeff_q2_err_CGM2"] = float(np.median(rel["S_CGM2"][m.q >= 30] * m.q[m.q >= 30] ** 2))
    s["thick_wall_coeff_q_nopt_minus_1"] = float(np.median((thick.n_CGM2 - 1) * thick.q))
    s["thick_wall_err_ratio_CGM2_over_CGM"] = float(np.median(rel["S_CGM2"][m.q >= 30] / rel["S_CGM"][m.q >= 30]))
    s["threshold_q_rel_error"] = {c: {lv: fnum(crossing(m.q, rel[c], v)) for lv, v in (("1%", .01), ("10%", .1), ("100%", 1.))}
                                  for c in rel}
    s["crossover_q_worse_than_FL_rel"] = {c: fnum(crossing(m.q, rel[c] / rel["S_FL"], 1.0)) for c in ("S_curv", "S_CGM", "S_CGM2")}
    s["crossover_q_worse_than_FL_log"] = {c: fnum(crossing(m.q, log[c] / log["S_FL"], 1.0)) for c in ("S_curv", "S_CGM", "S_CGM2")}
    fin = lambda c: m.q[np.isfinite(m[c]) & (m[c] > 0)]
    s["finiteness"] = {
        "S_CGM": {"predicted_threshold_q": 13 / 72, "smallest_q_finite": float(fin("S_CGM").min()),
                  "largest_q_infinite": fnum(m.q[~np.isfinite(m.S_CGM)].max()) if (~np.isfinite(m.S_CGM)).any() else None},
        "S_CGM2": {"predicted_threshold_q": 3 / 32, "smallest_q_finite": float(fin("S_CGM2").min()),
                   "largest_q_infinite": fnum(m.q[~np.isfinite(m.S_CGM2)].max()) if (~np.isfinite(m.S_CGM2)).any() else None},
    }
    s["S_CGM_exact_formula_max_rel_dev"] = float(((m.S_CGM - m.exact_CGM) / m.exact_CGM)[np.isfinite(m.exact_CGM)].abs().max())
    thin = m[(m.q <= 0.15) & m.S_num.notna()].sort_values("q")
    lq, lS = np.log(thin.q.values), np.log(thin.S_num.values)
    s["thin_wall_local_exponent_dlnS_dlnq"] = {f"{thin.q.values[i]:g}-{thin.q.values[i + 1]:g}":
                                               float((lS[i + 1] - lS[i]) / (lq[i + 1] - lq[i])) for i in range(len(thin) - 1)}
    s["thin_wall_ratio_S_num_over_S_TW"] = {f"{q:g}": float(S / (TW_COEFF / (l2 * q ** 3)))
                                            for q, S, l2 in zip(thin.q, thin.S_num, thin.l2)}
    s["M4_vs_M2_max_rel_diff"] = float(((m.S_M4 - m.S_M2).abs() / m.S_M2).max())
    s["M4_N_used"] = {f"{q:g}": int(n) for q, n in zip(m.q, m.N_M4)}
    s["table"] = [{"q": float(r.q), "S_num": float(r.S_num), "n_CGM2": fnum(r.n_CGM2),
                   **{f"rel_err_{c}": float(rel[c][i]) for c in rel}, **{f"log_err_{c}": float(log[c][i]) for c in log}}
                  for i, r in zip(m.index, m.itertuples())]
    return s


def sm_summary(sm, F_model):
    lines = {}
    for (presc, Mh), g in sm.groupby(["prescription", "Mh"]):
        g = g.sort_values("delta", ascending=False)
        ok = g[g.S_num.notna()].copy()
        ok["q"] = ok.lambda_min.abs() / ok.lambda2
        rel, log = errors(ok)
        d = {"points": int(len(g)), "solved": int(len(ok)), "status_counts": g.status.value_counts().to_dict(),
             "methods": g.method.value_counts().to_dict(),
             "planck_sensitive": int(g.planck_sensitive.astype(bool).sum()),
             "reduced_accuracy": int(g.reduced_accuracy.astype(bool).sum()),
             "log10_mu_min_range": [float(g.log10_mu_min.min()), float(g.log10_mu_min.max())],
             "unsolved_delta_GeV": [float(x) for x in g.delta[g.S_num.isna()]],
             "unsolved_q": [float(x) for x in (g.lambda_min.abs() / g.lambda2)[g.S_num.isna()]]}
        if len(ok) >= 4:
            near = ok[ok.delta <= 0.05]
            d["q_over_delta_per_GeV_near_critical"] = float(np.median(near.q / near.delta)) if len(near) else None
            d["delta_threshold_GeV"] = {c: {lv: fnum(crossing(ok.delta, rel[c], v)) for lv, v in (("1%", .01), ("10%", .1))}
                                        for c in rel}
            d["delta_crossover_worse_than_FL_rel"] = {c: fnum(crossing(ok.delta, rel[c] / rel["S_FL"], 1.0))
                                                      for c in ("S_curv", "S_CGM", "S_CGM2")}
            d["max_rel_error"] = {c: float(rel[c].max()) for c in rel}
            d["max_rel_error_finite"] = {c: float(rel[c][np.isfinite(rel[c])].max()) for c in rel}
            d["universality_FL_error_median_rel_dev_from_model"] = float(
                np.median(np.abs((ok.S_num / ok.S_FL - 1) / F_model(ok.q) - 1)))
            lx, ly = np.log(ok.delta.values), np.log(ok.S_num.values)
            d["dlnS_dlndelta_large_delta"] = float(np.polyfit(lx[:4], ly[:4], 1)[0])
            d["dlnS_dlndelta_small_delta"] = float(np.polyfit(lx[-4:], ly[-4:], 1)[0])
            d["table"] = [{"delta": float(r.delta), "q": float(r.q), "S_num": float(r.S_num), "status": r.status,
                           "planck_sensitive": bool(r.planck_sensitive), "reduced_accuracy": bool(r.reduced_accuracy),
                           **{f"rel_err_{c}": float(rel[c][i]) for c in rel}}
                          for i, r in zip(ok.index, ok.itertuples())]
        ps = g[g.planck_sensitive.astype(bool) | (g.status == "no-finite-bounce")]
        d["planck_affected_max_delta_GeV"] = float(ps.delta.max()) if len(ps) else None
        lines[f"{presc} Mh={Mh:g}"] = d
    return lines


def main():
    m = pd.read_csv(os.path.join(RES, "near_critical_model.csv")).sort_values("q", ascending=False).reset_index(drop=True)
    sm = pd.read_csv(os.path.join(RES, "near_critical_sm.csv"))
    s = {"model": model_summary(m)}
    mq, mF = np.log(m.q.values[::-1]), np.log((m.S_num / m.S_FL - 1).values[::-1])

    def F_model(q):          # model curve S/S_FL - 1 at the same q (log-log interpolation)
        return np.exp(np.interp(np.log(np.asarray(q, float)), mq, mF))
    s["sm_lines"] = sm_summary(sm, F_model)
    for name in ("controlled_cw1", "controlled_tree"):
        p = os.path.join(ROOT, "results", name + ".csv")
        p = p if os.path.exists(p) else p + ".gz"
        if os.path.exists(p):
            df = pd.read_csv(p, comment="#")
            d = df[df.S_num.notna() & (df.lambda2 > 0) & (df.planck_sensitive == 0) & (df.reduced_accuracy == 0)]
            qd = d.lambda_min.abs() / d.lambda2
            inside = (qd >= m.q.min()) & (qd <= m.q.max())
            dev = np.abs((d.S_num / d.S_FL - 1)[inside] / F_model(qd[inside]) - 1)
            s[f"universality_{name}"] = {"n": int(inside.sum()), "median_rel_dev": float(np.median(dev)),
                                         "p90_rel_dev": float(np.quantile(dev, 0.9)), "max_rel_dev": float(dev.max()),
                                         "q_range": [float(qd.min()), float(qd.max())]}
    with open(os.path.join(RES, "near_critical_summary.json"), "w") as fh:
        json.dump(s, fh, indent=1, default=float)
    short = {"model": {k: v for k, v in s["model"].items() if k != "table"},
             "sm_lines": {k: {kk: vv for kk, vv in v.items() if kk != "table"} for k, v in s["sm_lines"].items()},
             **{k: v for k, v in s.items() if k.startswith("universality")}}
    print(json.dumps(short, indent=1, default=float))


if __name__ == "__main__":
    main()
