"""Generate every figure in figures/ from the CSV outputs (no data are modified before plotting).

Inputs: results/controlled_cw1.csv, results/controlled_tree.csv, results/exploratory_cw1.csv,
        validation/results/method_comparison.csv, validation/results/convergence.csv,
        validation/results/legacy_comparison.csv (optional)
Usage:  python validation/make_figures.py
"""
from __future__ import annotations

import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11, "legend.fontsize": 8.5,
                     "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight"})

STATUS = ["stable", "metastable", "unstable", "non-perturbative", "no-bounce", "solver-failure", "invalid-input"]
STATUS_COLORS = ["#1b9e77", "#fdb863", "#b2182b", "#636363", "#cfcfcf", "#000000", "#ffffff"]
REF_POINT = (125.15, 173.34)   # central inputs of Buttazzo et al. (2013)


def read(path):
    p = os.path.join(ROOT, path)
    if not os.path.exists(p):
        return None, {}
    meta = {}
    with open(p) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            k, _, v = line[1:].partition("=")
            meta[k.strip()] = v.strip()
    return pd.read_csv(p, comment="#"), meta


def grid(df, col):
    mh = np.sort(df.Mh.unique())
    mt = np.sort(df.Mt.unique())
    Z = df.pivot_table(index="Mt", columns="Mh", values=col, aggfunc="first").reindex(index=mt, columns=mh).values
    return mh, mt, Z


def status_map(ax, df, title, legend=True):
    code = df.status.map({s: i for i, s in enumerate(STATUS)})
    d = df.assign(code=code)
    mh, mt, Z = grid(d, "code")
    cmap = ListedColormap(STATUS_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, len(STATUS) + 0.5), cmap.N)
    ax.pcolormesh(mh, mt, Z, cmap=cmap, norm=norm, shading="nearest", rasterized=True)
    if legend:
        present = [s for s in STATUS if (df.status == s).any()]
        handles = [plt.Rectangle((0, 0), 1, 1, fc=STATUS_COLORS[STATUS.index(s)], ec="k", lw=0.3) for s in present]
        ax.legend(handles, present, loc="lower right", frameon=True, title="numerical reference")
    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    ax.set_title(title)


def fig_phase_diagram(df, meta):
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    status_map(ax, df, "SM vacuum stability (controlled region)", legend=False)
    mh, mt, Lm = grid(df, "lambda_min")
    _, _, Pn = grid(df, "log10_p_num")
    _, _, Pl = grid(df, "log10_p_FL_lit")
    _, _, Pb = grid(df, "log10_p_FL_bscale")
    c1 = ax.contour(mh, mt, Lm, levels=[0.0], colors="k", linewidths=1.6)
    c2 = ax.contour(mh, mt, Pn, levels=[0.0], colors="k", linewidths=1.6, linestyles="-")
    c3 = ax.contour(mh, mt, Pl, levels=[0.0], colors="#2166ac", linewidths=1.3, linestyles="--")
    c4 = ax.contour(mh, mt, Pb, levels=[0.0], colors="#2166ac", linewidths=1.0, linestyles=":")
    ps = df[df.planck_sensitive == 1]
    if len(ps):
        ax.scatter(ps.Mh, ps.Mt, s=2, c="k", marker=".", label="Planck-sensitive", rasterized=True)
    ax.plot(*REF_POINT, marker="*", ms=12, mfc="w", mec="k", ls="none")
    ax.annotate("Buttazzo et al. (2013)\ncentral inputs", REF_POINT, xytext=(6, -14), textcoords="offset points", fontsize=8)
    from matplotlib.lines import Line2D
    handles = [plt.Rectangle((0, 0), 1, 1, fc=STATUS_COLORS[i], ec="k", lw=0.3) for i in (0, 1, 2)]
    handles += [Line2D([], [], color="k", lw=1.6), Line2D([], [], color="#2166ac", ls="--"),
                Line2D([], [], color="#2166ac", ls=":")]
    labels = ["stable", "metastable", "unstable", r"boundaries: $\lambda_{\min}=0$ and $\wp=1$ (numerical bounce)",
              r"$\wp=1$: $S_{\rm FL}$, $\Lambda_B=\arg\min\lambda$ (literature)", r"$\wp=1$: $S_{\rm FL}$, $\Lambda_B=1/R$"]
    ax.legend(handles, labels, loc="upper left", fontsize=7.5, frameon=True)
    ax.text(0.99, 0.01, f"prescription {meta.get('prescription','?')}, matching {meta.get('matching','?')},\n"
            rf"$\alpha_s(M_Z)$={meta.get('alpha3','?')}, $H_0$={meta.get('H0_km_s_Mpc','?')} km/s/Mpc",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, bbox=dict(fc="w", ec="0.7", lw=0.5))
    fig.savefig(os.path.join(FIG, "phase_diagram_controlled.png"))
    plt.close(fig)


def fig_discrepancy_maps(df):
    d = df[df.S_num.notna()].copy()
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.3), sharey=True)
    specs = [("S_FL", r"$(S_{\rm num}-S_{\rm FL})/S_{\rm num}$ (literature; rigorous lower bound)"),
             ("S_curv", r"$|S_{\rm curv}-S_{\rm num}|/S_{\rm num}$ (curvature-corrected)"),
             ("S_CGM", r"$(S_{\rm CGM}-S_{\rm num})/S_{\rm num}$ (variational upper estimate)")]
    norm = LogNorm(1e-9, 1e0)
    for ax, (col, title) in zip(axs, specs):
        d["err"] = (d[col] - d.S_num).abs() / d.S_num
        mh, mt, Z = grid(df.merge(d[["Mh", "Mt", "err"]], on=["Mh", "Mt"], how="left"), "err")
        m = ax.pcolormesh(mh, mt, np.clip(Z, 1e-9, 1), norm=norm, cmap="viridis", shading="nearest", rasterized=True)
        cs = ax.contour(mh, mt, Z, levels=[1e-6, 1e-4, 1e-3, 1e-2, 1e-1], colors="w", linewidths=0.7)
        ax.clabel(cs, fmt=lambda v: f"{v:.0e}", fontsize=7)
        _, _, Lm = grid(df, "lambda_min")
        ax.contour(mh, mt, Lm, levels=[0.0], colors="r", linewidths=1.2)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(r"$M_h$ [GeV]")
    axs[0].set_ylabel(r"$M_t$ [GeV]")
    fig.colorbar(m, ax=axs, shrink=0.9, label="relative deviation from the numerical reference")
    fig.text(0.01, 0.01, r"Red line: absolute-stability boundary $\lambda_{\min}=0$. White: no bounce (stable).", fontsize=8)
    fig.savefig(os.path.join(FIG, "analytic_vs_numerical_maps.png"))
    plt.close(fig)


def fig_error_scaling(dfs, comp):
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    styles = {"cw1": ("o", 0.25), "tree": ("s", 0.25)}
    for name, df in dfs.items():
        d = df[(df.S_num.notna()) & (df.lambda2 > 0) & (df.planck_sensitive == 0)]
        q = d.lambda_min.abs() / d.lambda2
        mk, al = styles.get(name, ("o", 0.3))
        ax.loglog(q, (d.S_num - d.S_FL) / d.S_num, mk, ms=2, alpha=al, color="#b2182b", rasterized=True)
        ax.loglog(q, (d.S_curv - d.S_num).abs() / d.S_num, mk, ms=2, alpha=al, color="#2166ac", rasterized=True)
        ax.loglog(q, (d.S_CGM - d.S_num).abs() / d.S_num, mk, ms=2, alpha=al, color="#1b9e77", rasterized=True)
    if comp is not None:
        c = comp[comp.case.str.startswith("B2/3")]
        ref = c[c.method == "M4"].set_index("case")
        for case, row in ref.iterrows():
            ax.loglog(row.q, (row.S - row.S_FL) / row.S, "k*", ms=9, zorder=5)
    qq = np.logspace(-1, 3.5, 200)
    ax.loglog(qq, (13 / 72) / qq, "k--", lw=1, label=r"first order: $\frac{13}{72}\,\lambda''/|\lambda_{\min}|$")
    ax.plot([], [], "o", color="#b2182b", label=r"$S_{\rm FL}$ (literature)")
    ax.plot([], [], "o", color="#2166ac", label=r"$S_{\rm curv}$")
    ax.plot([], [], "o", color="#1b9e77", label=r"$S_{\rm CGM}$")
    ax.plot([], [], "k*", ms=8, label="parabolic model (benchmark)")
    ax.set_xlabel(r"$q=|\lambda_{\min}|/\lambda''_{\min}$  (thin wall $\leftarrow$ $\rightarrow$ thick wall)")
    ax.set_ylabel("relative deviation from numerical reference")
    ax.set_ylim(1e-12, 3)
    ax.legend(loc="lower left")
    ax.set_title("Accuracy of analytical approximations vs. conformal-breaking parameter")
    ax.grid(alpha=0.3, which="both", lw=0.3)
    fig.savefig(os.path.join(FIG, "analytic_error_scaling.png"))
    plt.close(fig)


def fig_method_comparison(comp):
    cases = [c for c in comp.case.unique() if not c.startswith("B1")]
    methods = ["M0", "M1", "M2", "M3", "M4"]
    colors = {"M0": "#999999", "M1": "#e66101", "M2": "#5e3c99", "M3": "#fdb863", "M4": "#1b9e77"}
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    x = np.arange(len(cases))
    for j, m in enumerate(methods):
        dev, tim = [], []
        for c in cases:
            sub = comp[comp.case == c]
            good = sub[(sub.method.isin(["M1", "M2", "M4"])) & (sub.ok == True)]  # noqa: E712
            ref = np.median(good.S) if len(good) else np.nan
            r = sub[sub.method == m]
            ok = len(r) and bool(r.ok.iloc[0]) and np.isfinite(r.S.iloc[0]) and r.S.iloc[0] < 1e99
            dev.append(abs(r.S.iloc[0] - ref) / ref if ok else np.nan)
            tim.append(r.seconds.iloc[0] if len(r) else np.nan)
            if not ok:
                a1.plot(x[len(dev) - 1] + (j - 2) * 0.15, 3e0, "x", color=colors[m], ms=7)
        a1.semilogy(x + (j - 2) * 0.15, np.maximum(dev, 1e-16), "o", color=colors[m], label=m)
        a2.semilogy(x + (j - 2) * 0.15, tim, "s", color=colors[m])
    a1.axhline(1e-10, color="k", lw=0.5, ls=":")
    a1.set_ylabel("|S - consensus| / consensus\n(x = failure / sentinel / negative)")
    a1.set_ylim(1e-16, 1e1)
    a1.legend(ncol=5, loc="upper center", fontsize=8)
    a2.set_ylabel("wall time [s] (Python)")
    a2.set_xticks(x)
    a2.set_xticklabels([c.replace("B2/3 ", "").replace("B4 SM ", "SM ").replace("B5 SM ", "SM* ") for c in cases],
                       rotation=55, ha="right", fontsize=7.5)
    a1.set_title("Candidate bounce solvers on the benchmark set (consensus of M1, M2, M4)")
    fig.savefig(os.path.join(FIG, "method_comparison.png"))
    plt.close(fig)


def fig_convergence(conv):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for case, sub in conv.groupby("case", sort=False):
        ref = np.median(sub[(sub.parameter == "rtol") & (sub.value.astype(str) == "1e-13") & (sub.ok == True)].S)  # noqa: E712
        m4 = sub[(sub.method == "M4") & (sub.parameter == "N")].copy()
        m4["N"] = m4.value.astype(float)
        a1.semilogy(m4.N, np.maximum((m4.S - ref).abs() / ref, 1e-17), "o-", ms=3, label=case)
        for m, ls in (("M2", "-"), ("M1", "--")):
            s = sub[(sub.method == m) & (sub.parameter == "rtol")].copy()
            s["rt"] = s.value.astype(float)
            a2.loglog(s.rt, np.maximum((s.S - ref).abs() / ref, 1e-17), "o" + ls, ms=3, label=f"{m}: {case}" if m == "M2" else None)
    a1.set_xlabel("Chebyshev nodes N (M4)")
    a1.set_ylabel("|S - reference| / reference")
    a1.set_title("Spectral collocation: convergence in N")
    a1.legend(fontsize=7)
    a2.set_xlabel("integrator tolerance rtol (M2 solid, M1 dashed)")
    a2.set_title("Shooting: convergence in tolerance")
    a2.legend(fontsize=7)
    for a in (a1, a2):
        a.grid(alpha=0.3, which="both", lw=0.3)
    fig.savefig(os.path.join(FIG, "convergence.png"))
    plt.close(fig)


def fig_profiles():
    sys.path.insert(0, os.path.join(ROOT, "python"))
    from smvac.reference import analytic as A, bounce as B, potential as P, sm
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    cases = [(173.34, r"$M_t=173.34$ GeV: thick wall"), (171.119296 + 0.01, r"$M_t=M_t^{\rm crit}+0.01$ GeV: thin wall")]
    for ax, (Mt, label) in zip(axs, cases):
        pot = P.sm_potential(sm.Inputs(125.15, Mt), "cw1")
        an = A.estimates(pot)
        r = None
        for N in (256, 512, 1024):
            r = B.solve_m4(pot, an.lambda_min, an.lnh_CGM, N=N)
            if r.ok and abs(r.derrick) < 1e-9 and abs(r.dilatation) < 1e-9:
                break
        tau, f = r.profile
        rr, phi = np.exp(tau), f * np.exp(-tau)
        phi0 = math.exp(r.lnphi0)
        Rfl = math.sqrt(8.0 / abs(an.lambda_min)) / phi0      # FL profile with lambda_min and the same phi(0)
        ax.loglog(rr, phi, "-", color="#1b9e77", lw=2, label="numerical bounce (M4)")
        ax.loglog(rr, phi0 * Rfl ** 2 / (rr ** 2 + Rfl ** 2), "--", color="0.3",
                  label=r"Fubini-Lipatov, $\lambda_{\min}$, same $\phi(0)$")
        ax.axhline(math.exp(an.x_min), color="#b2182b", lw=0.8, ls=":", label=r"field at $\lambda_{\min}$")
        ax.set_title(label)
        ax.set_xlabel(r"$r$ [GeV$^{-1}$]")
        ax.set_ylabel(r"$\phi(r)$ [GeV]")
        ax.set_ylim(phi0 * 1e-6, phi0 * 3)
        ax.text(0.03, 0.05, rf"$S_{{\rm num}}$ = {r.S:.7g}" + "\n" + rf"$S_{{\rm FL}}$ = {an.S_FL:.7g}",
                transform=ax.transAxes, fontsize=8, bbox=dict(fc="w", ec="0.7", lw=0.5))
        ax.legend(fontsize=7.5, loc="upper right")
        ax.grid(alpha=0.3, which="both", lw=0.3)
    fig.suptitle(r"Bounce profiles at $M_h=125.15$ GeV (prescription cw1)", fontsize=11)
    fig.savefig(os.path.join(FIG, "bounce_profiles.png"))
    plt.close(fig)


def fig_exploratory(df):
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    status_map(ax, df, "EXPLORATORY: uncontrolled extrapolation of the inputs")
    ax.add_patch(plt.Rectangle((120, 166), 10, 16, fill=False, ec="b", lw=1.5))
    ax.text(131, 160, "controlled\nregion", color="b", fontsize=8)
    ax.text(0.5, 0.55, "NOT A PHYSICAL PREDICTION", transform=ax.transAxes, ha="center", fontsize=14, color="k",
            alpha=0.25, rotation=30)
    fig.savefig(os.path.join(FIG, "exploratory_wide_scan.png"))
    plt.close(fig)


def main():
    cw1, meta = read("results/controlled_cw1.csv")
    tree, _ = read("results/controlled_tree.csv")
    comp, _ = read("validation/results/method_comparison.csv")
    conv, _ = read("validation/results/convergence.csv")
    expl, _ = read("results/exploratory_cw1.csv")
    if cw1 is not None:
        fig_phase_diagram(cw1, meta)
        fig_discrepancy_maps(cw1)
        fig_error_scaling({k: v for k, v in (("cw1", cw1), ("tree", tree)) if v is not None}, comp)
    if comp is not None:
        fig_method_comparison(comp)
    if conv is not None:
        fig_convergence(conv)
    if expl is not None:
        fig_exploratory(expl)
    fig_profiles()
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
