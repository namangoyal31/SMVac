"""Print Markdown tables summarising validation/results/method_comparison.csv (used in
docs/NUMERICAL_METHOD_COMPARISON.md). The consensus value per case is the median of the successful runs of
M1, M2 and M4 (three structurally different methods)."""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.stdout.reconfigure(encoding="utf-8")
OUT = open(os.path.join(ROOT, "validation", "results", "method_comparison_tables.md"), "w", encoding="utf-8")
_print = print


def print(*a, **k):  # echo to stdout and to the Markdown file used by docs/NUMERICAL_METHOD_COMPARISON.md
    _print(*a, **k)
    _print(*a, file=OUT)

df = pd.read_csv(os.path.join(ROOT, "validation", "results", "method_comparison.csv"))
methods = ["M0", "M1", "M2", "M3", "M4"]


def fmt(x):
    if not np.isfinite(x):
        return "fail"
    return "0" if x == 0 else f"{x:.0e}"


rows_dev, rows_an = [], []
for case, sub in df.groupby("case", sort=False):
    good = sub[sub.method.isin(["M1", "M2", "M4"]) & (sub.ok == True)]  # noqa: E712
    ref = float(np.median(good.S)) if len(good) else np.nan
    spread = float((good.S.max() - good.S.min()) / ref) if len(good) > 1 else np.nan
    line = [case, f"{ref:.12g}", fmt(spread)]
    for m in methods:
        r = sub[sub.method == m]
        if not len(r) or not bool(r.ok.iloc[0]) or not np.isfinite(r.S.iloc[0]) or r.S.iloc[0] >= 1e99:
            line.append("fail" if len(r) else "–")
        elif r.S.iloc[0] < 0:
            line.append(f"{r.S.iloc[0]:.1e} (negative)")
        else:
            d = abs(r.S.iloc[0] - ref) / ref
            line.append(f"{fmt(d)} / {r.seconds.iloc[0]:.1f}s")
    rows_dev.append(line)
    a = sub.iloc[0]
    q = a.q if np.isfinite(a.q) else np.nan
    rows_an.append([case, f"{q:.3g}" if np.isfinite(q) else "–", f"{(ref - a.S_FL) / ref:+.2e}",
                    f"{(a.S_curv - ref) / ref:+.2e}" if np.isfinite(a.S_curv) else "–",
                    f"{(a.S_CGM - ref) / ref:+.2e}" if np.isfinite(a.S_CGM) else "–"])

print("| case | consensus S | spread (M1,M2,M4) | " + " | ".join(f"{m}: dev / time" for m in methods) + " |")
print("|" + "---|" * (3 + len(methods)))
for r in rows_dev:
    print("| " + " | ".join(r) + " |")
print()
print("| case | q = λ_min/λ″ (abs.) | (S−S_FL)/S | (S_curv−S)/S | (S_CGM−S)/S |")
print("|---|---|---|---|---|")
for r in rows_an:
    print("| " + " | ".join(r) + " |")
