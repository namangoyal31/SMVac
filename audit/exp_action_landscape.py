"""E4: landscape of the legacy objective S(R) minimised by golden-section search.

For benchmark points, evaluates the legacy action_at_R on a fine grid of ln R across the
legacy bracket [ln R_opt - 8, ln R_opt + 8] and records kinetic/potential parts, the legacy
GSS result and the analytic estimate S_approx = 8 pi^2 / (3 |lambda_eff,min|).
Outputs: audit/results/action_landscape.csv and audit/results/action_landscape.png
"""
import math, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from audit import legacy_port as L

POINTS = [(125.1, 173.1), (115.0, 180.0), (134.75, 176.5), (5.0, 105.0), (5.0, 96.0)]
rows = []
fig, axes = plt.subplots(1, len(POINTS), figsize=(4.2 * len(POINTS), 3.8))
for ax, (Mh, Mt) in zip(axes, POINTS):
    res = L.classify(Mh, Mt)
    run = res['run']
    mu_opt = math.exp(run.t_min_lambda / 2)
    lr0 = math.log(run.mu_inst / mu_opt)
    grid = np.linspace(lr0 - 8, lr0 + 8, 801)
    S, T, U = [], [], []
    for lr in grid:
        v = L.action_at_R(run.table, run.mu_inst, math.exp(lr), return_parts=True)
        if isinstance(v, tuple):
            S.append(v[0]); T.append(v[1]); U.append(v[2])
        else:
            S.append(v); T.append(np.nan); U.append(np.nan)
        rows.append(dict(Mh=Mh, Mt=Mt, lnR=lr, log10_muR=math.log10(run.mu_inst) - lr / math.log(10),
                         S=S[-1], kinetic=T[-1], potential=U[-1]))
    S = np.array(S)
    muR = np.log10(run.mu_inst) - grid / np.log(10)
    finite = S < 1e99
    ax.plot(muR[finite], S[finite], lw=1, label='legacy $S(R)$')
    ax.axhline(res['S_approx'], color='k', ls='--', lw=0.8, label=r'$8\pi^2/3|\lambda_{\rm eff,min}|$')
    ax.plot([math.log10(res['mu_R_min'])], [res['S_exact']], 'r*', ms=10, label='legacy GSS result')
    ax.set_title(f"$M_h$={Mh}, $M_t$={Mt}\nS_exact={res['S_exact']:.6g}", fontsize=9)
    ax.set_xlabel(r'$\log_{10}(\mu_R/{\rm GeV})$, $\mu_R=1/R$')
    lo, hi = np.nanpercentile(S[finite], [1, 99]) if finite.any() else (0, 1)
    ax.set_ylim(min(lo, res['S_exact'], res['S_approx']) - 0.1 * abs(hi - lo), max(hi, res['S_approx']) + 0.1 * abs(hi - lo))
    ax.grid(alpha=0.3)
    n_sentinel = int((~finite).sum()); n_neg = int((S[finite] < 0).sum())
    ax.text(0.02, 0.02, f"sentinel 1e100: {n_sentinel}/801\nS<0: {n_neg}/801", transform=ax.transAxes, fontsize=7)
axes[0].set_ylabel('legacy action functional')
axes[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig('audit/results/action_landscape.png', dpi=150)
with open('audit/results/action_landscape.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote audit/results/action_landscape.{csv,png}")
