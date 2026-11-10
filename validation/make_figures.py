"""Generate every figure in figures/ from the CSV outputs (no data are modified before plotting).

Inputs: results/controlled_cw1.csv, results/controlled_cw1_fallback.csv, results/controlled_tree.csv,
        results/exploratory_cw1.csv, validation/results/method_comparison.csv, validation/results/convergence.csv,
        validation/results/legacy_comparison.csv, validation/results/near_critical_model.csv,
        validation/results/near_critical_sm.csv (CSV files may be gzip-compressed, *.csv.gz).
Two figures also solve bounces with the Python reference (bounce_profiles.png, variational_failure.png).
Captions: figures/README.md.
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
STATUS_LABELS = {
    "stable": r"absolutely stable ($\lambda_{\min}\geq0$)",
    "metastable": r"metastable ($\wp<1$)",
    "unstable": r"unstable ($\wp>1$)",
    "non-perturbative": r"non-perturbative ($|\lambda|$ or $y_t>4\pi$ below $M_{\rm Pl}$)",
    "no-bounce": r"no finite bounce: $\lambda_{\min}$ at or just below $M_{\rm Pl}$ (Planck plateau)",
    "solver-failure": "unresolved (solver failure)",
    "invalid-input": "invalid input (matching outside its domain)",
}


def read(path):
    p = os.path.join(ROOT, path)
    if not os.path.exists(p) and os.path.exists(p + ".gz"):
        p += ".gz"                     # committed scan outputs are gzip-compressed
    if not os.path.exists(p):
        return None, {}
    meta = {}
    import gzip
    with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as fh:
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
        ax.legend(handles, [STATUS_LABELS[p] for p in present], loc="upper center", bbox_to_anchor=(0.5, -0.12),
                  ncol=2, frameon=True, fontsize=7, title="numerical reference")
    ax.set_xlabel(r"$M_h$ [GeV]")
    ax.set_ylabel(r"$M_t$ [GeV]")
    ax.set_title(title)


def apply_fallback(df, fb):
    """Replace solver-failure rows by the accepted shooting-fallback results (flagged in column 'fallback')."""
    df = df.copy()
    df["fallback"] = 0
    if fb is None:
        return df
    for r in fb.itertuples():
        i = np.nonzero((np.isclose(df.Mh, r.Mh)) & (np.isclose(df.Mt, r.Mt)))[0]
        if not len(i):
            continue
        j = i[0]
        if r.fallback_status in ("metastable", "unstable"):
            df.loc[j, ["status", "S_num", "log10_p_num", "planck_sensitive", "fallback"]] = [
                r.fallback_status, r.S_fallback, r.log10_p, int(r.planck_sensitive), 1]
        elif r.fallback_status == "no-finite-bounce":   # Planck plateau: same physics as the no-bounce points
            df.loc[j, ["status", "planck_sensitive", "fallback"]] = ["no-bounce", 1, 1]
    return df


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
    shown = [i for i, s in enumerate(STATUS) if (df.status == s).any()]
    handles = [plt.Rectangle((0, 0), 1, 1, fc=STATUS_COLORS[i], ec="k", lw=0.3) for i in shown]
    labels = [STATUS_LABELS[STATUS[i]] for i in shown]
    if len(ps):
        handles.append(Line2D([], [], color="k", marker=".", ls="none"))
        labels.append(r"Planck-sensitive ($\phi(0)>M_{\rm Pl}$ or $\lambda_{\min}$ at/near $M_{\rm Pl}$)")
    handles += [Line2D([], [], color="k", lw=1.6), Line2D([], [], color="#2166ac", ls="--"),
                Line2D([], [], color="#2166ac", ls=":")]
    labels += [r"boundaries: $\lambda_{\min}=0$ and $\wp=1$ (numerical bounce)",
               r"$\wp=1$: $S_{\rm FL}$, $\Lambda_B=\arg\min\lambda$ (literature)", r"$\wp=1$: $S_{\rm FL}$, $\Lambda_B=1/R$"]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=7, frameon=True)
    ax.text(0.99, 0.01, f"prescription {meta.get('prescription','?')}, matching {meta.get('matching','?')},\n"
            rf"$\alpha_s(M_Z)$={meta.get('alpha3','?')}, $H_0$={meta.get('H0_km_s_Mpc','?')} km/s/Mpc",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, bbox=dict(fc="w", ec="0.7", lw=0.5))
    fig.savefig(os.path.join(FIG, "phase_diagram_controlled.png"))
    plt.close(fig)


def fig_discrepancy_maps(df):
    d = df[df.S_num.notna()].copy()
    d["S_E1"] = d.S_FL * (1 + (13 / 72) * d.lambda2 / d.lambda_min.abs())   # Espinosa (2020) eq. (3.21)
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.3), sharey=True)
    specs = [("S_FL", r"$(S_{\rm num}-S_{\rm FL})/S_{\rm num}$ (FL; rigorous lower bound)"),
             ("S_E1", r"$|S_{\rm E1}-S_{\rm num}|/S_{\rm num}$ (first order, Espinosa 2020)"),
             ("S_curv", r"$|S_{\rm curv}-S_{\rm num}|/S_{\rm num}$ (resummed first order)"),
             ("S_CGM", r"$(S_{\rm CGM}-S_{\rm num})/S_{\rm num}$ (variational, FL trial)")]
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
    fig.text(0.01, -0.02, r"Red line: absolute-stability boundary $\lambda_{\min}=0$. White: no finite bounce "
             r"(absolutely stable, or Planck plateau at $M_h\gtrsim127.5$ GeV). Prescription cw1; the 78 Planck-sensitive "
             "shooting-fallback points are included.", fontsize=8)
    fig.savefig(os.path.join(FIG, "analytic_vs_numerical_maps.png"))
    plt.close(fig)


def fig_error_scaling(dfs, comp):
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    styles = {"cw1": ("o", 0.25), "tree": ("s", 0.25)}
    for name, df in dfs.items():
        d = df[(df.S_num.notna()) & (df.lambda2 > 0) & (df.planck_sensitive == 0) & (df.reduced_accuracy == 0)]
        q = d.lambda_min.abs() / d.lambda2
        mk, al = styles.get(name, ("o", 0.3))
        ax.loglog(q, (d.S_num - d.S_FL) / d.S_num, mk, ms=2, alpha=al, color="#b2182b", rasterized=True)
        ax.loglog(q, (d.S_FL * (1 + (13 / 72) / q) - d.S_num).abs() / d.S_num, mk, ms=2, alpha=al, color="#e6ab02",
                  rasterized=True)
        ax.loglog(q, (d.S_curv - d.S_num).abs() / d.S_num, mk, ms=2, alpha=al, color="#2166ac", rasterized=True)
        ax.loglog(q, (d.S_CGM - d.S_num).abs() / d.S_num, mk, ms=2, alpha=al, color="#1b9e77", rasterized=True)
    if comp is not None:
        c = comp[comp.case.str.startswith("B2/3")]
        ref = c[c.method == "M4"].set_index("case")
        for case, row in ref.iterrows():
            ax.loglog(row.q, (row.S - row.S_FL) / row.S, "k*", ms=9, zorder=5)
    qq = np.logspace(-1, 3.5, 200)
    ax.loglog(qq, (13 / 72) / qq, "k--", lw=1, label=r"thick-wall law $\frac{13}{72}q^{-1}$ (Espinosa 2020)")
    ax.plot([], [], "o", color="#b2182b", label=r"$S_{\rm FL}$ (literature; lower bound)")
    ax.plot([], [], "o", color="#e6ab02", label=r"$S_{\rm E1}=S_{\rm FL}(1+\frac{13}{72q})$ (first order)")
    ax.plot([], [], "o", color="#2166ac", label=r"$S_{\rm curv}=S_{\rm FL}/(1-\frac{13}{72q})$")
    ax.plot([], [], "o", color="#1b9e77", label=r"$S_{\rm CGM}$ (variational, FL trial)")
    ax.plot([], [], "k*", ms=8, label="parabolic model (benchmark)")
    ax.set_xlabel(r"$q=|\lambda_{\min}|/\lambda''_{\min}$  (thin wall $\leftarrow$ $\rightarrow$ thick wall)")
    ax.set_ylabel("relative deviation from numerical reference")
    ax.set_ylim(1e-8, 3)
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
            if not len(r):          # not run for this case (e.g. shooting on the degenerate constant case)
                dev.append(np.nan)
                tim.append(np.nan)
                continue
            S = r.S.iloc[0]
            ok = bool(r.ok.iloc[0]) and np.isfinite(S) and 0 < S < 1e99
            dev.append(abs(S - ref) / ref if ok else np.nan)
            tim.append(r.seconds.iloc[0])
            if not ok:   # failure, sentinel or negative "action"
                a1.plot(x[len(dev) - 1] + (j - 2) * 0.15, 3e0, "x", color=colors[m], ms=7, mew=2)
        a1.semilogy(x + (j - 2) * 0.15, np.maximum(dev, 1e-16), "o", color=colors[m], label=m)
        a2.semilogy(x + (j - 2) * 0.15, tim, "s", color=colors[m])
    a1.axhline(1e-10, color="k", lw=0.5, ls=":")
    a1.set_ylabel("|S - consensus| / consensus\n(x = failure / sentinel / negative)")
    a1.set_ylim(1e-16, 1e1)
    a1.legend(ncol=5, loc="center left", fontsize=8, framealpha=0.9)
    a1.text(0.995, 0.03, "points on the lower edge coincide with the consensus (median member)", transform=a1.transAxes,
            ha="right", fontsize=7, color="0.3")
    a2.set_ylabel("wall time [s]\n(Python, machine under load)")
    a2.set_xticks(x)
    a2.set_xticklabels([c.replace("B2/3 ", "").replace("B4 SM ", "SM ").replace("B5 SM ", "SM* ") for c in cases],
                       rotation=55, ha="right", fontsize=7.5)
    a1.set_title("Candidate bounce solvers on the benchmark set (consensus of M1, M2, M4)")
    fig.savefig(os.path.join(FIG, "method_comparison.png"))
    plt.close(fig)


def fig_convergence(conv):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for k, (case, sub) in enumerate(conv.groupby("case", sort=False)):
        c = colors[k % len(colors)]
        ref = np.median(sub[(sub.parameter == "rtol") & (sub.value.astype(str) == "1e-13") & (sub.ok == True)].S)  # noqa: E712
        m4 = sub[(sub.method == "M4") & (sub.parameter == "N")].copy()
        m4["N"] = m4.value.astype(float)
        a1.semilogy(m4.N, np.maximum((m4.S - ref).abs() / ref, 1e-17), "o-", ms=3, color=c, label=case)
        for m, ls in (("M2", "-"), ("M1", "--")):
            s = sub[(sub.method == m) & (sub.parameter == "rtol")].copy()
            s["rt"] = s.value.astype(float)
            a2.loglog(s.rt, np.maximum((s.S - ref).abs() / ref, 1e-17), "o" + ls, ms=3, color=c,
                      label=case if m == "M2" else None)
    a1.text(0.98, 0.45, "plateau ~5e-12 = accuracy of the shooting\nreference, not of M4 (VALIDATION_REPORT 4.1)",
            transform=a1.transAxes, ha="right", fontsize=7, color="0.3")
    a1.set_xlabel("Chebyshev nodes N (M4)")
    a1.set_ylabel("|S - reference| / reference\n(reference: shooting consensus at rtol 1e-13)")
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
    cases = [(173.34, r"$M_t=173.34$ GeV (central inputs)"), (171.119296 + 0.01, r"$M_t=M_t^{\rm crit}+0.01$ GeV")]
    for ax, (Mt, label) in zip(axs, cases):
        pot = P.sm_potential(sm.Inputs(125.15, Mt), "cw1")
        an = A.estimates(pot)
        label += rf": $q={abs(an.lambda_min) / an.lambda2:.3g}$"
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
        ax.set_xlim(math.exp(r.lnR) * 10 ** -3, math.exp(r.lnR) * 10 ** 3)
        ax.text(0.03, 0.30, rf"$S_{{\rm num}}$ = {r.S:.7g}" + "\n" + rf"$S_{{\rm FL}}$ = {an.S_FL:.7g}",
                transform=ax.transAxes, fontsize=8, bbox=dict(fc="w", ec="0.7", lw=0.5))
        ax.legend(fontsize=7.5, loc="lower left")
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


def fig_legacy(legacy):
    """Why the legacy construction fails: (a) its objective on a smooth model potential, decomposed into the FL
    term 8pi^2/(3|lambda_c|) and the first-order correction ~ 1/lambda_c^2; (b) legacy vs new on the SM grid."""
    sys.path.insert(0, os.path.join(ROOT, "python"))
    from smvac.reference import analytic as A, bounce as B, potential as P
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    l2, q = 5e-5, 25.0
    pot = P.parabola_model(-q * l2, l2, 39.0)
    an = A.estimates(pot, (25, 55))
    ref = B.solve_m4(pot, an.lambda_min, an.lnh_CGM, N=384).S
    s_ = np.linspace(-30, 30, 6001)
    w = np.full(s_.size, 2.0); w[1::2] = 4.0; w[0] = w[-1] = 1.0; w *= (s_[1] - s_[0]) / 3
    lnR = np.linspace(-46.0, -31.95, 1500)      # 1/R across the whole window where lambda(1/R) < 0
    lead, corr, lamc = [], [], []
    for lr in lnR:
        lc = float(pot.lam(-lr))
        lamc.append(lc)
        if lc >= 0:
            lead.append(np.nan); corr.append(np.nan); continue
        Af = math.sqrt(2 / abs(lc))
        f = Af / np.cosh(s_)
        x = np.log(f) - (s_ + lr)
        lead.append(8 * math.pi ** 2 / (3 * abs(lc)))
        corr.append(0.5 * math.pi ** 2 * np.dot(w, (pot.lam(x) - lc) * f ** 4))
    lead, corr = np.array(lead), np.array(corr)
    a1.plot(-lnR / math.log(10), lead + corr, "k-", lw=1.8, label=r"legacy objective $S_{\rm legacy}(R)$")
    a1.plot(-lnR / math.log(10), lead, "--", color="#2166ac", label=r"FL term $8\pi^2/3|\lambda_c|$, $\lambda_c=\lambda(1/R)$")
    a1.plot(-lnR / math.log(10), corr, ":", color="#b2182b", label=r"first-order correction $\propto\delta\lambda/\lambda_c^2$")
    a1.axhline(ref, color="#1b9e77", lw=1.2, label=f"numerical bounce S = {ref:.6g}")
    tot = lead + corr
    i = int(np.nanargmin(tot))
    a1.plot(-lnR[i] / math.log(10), tot[i], "kv", ms=7)
    a1.annotate(r"$S_{\rm legacy}\to-\infty$ as $\lambda(1/R)\to0^-$" + f"\n({tot[i]:.2g} at the last grid point;\n"
                r"the legacy minimiser returns $-3.0\times10^{17}$)", (-lnR[i] / math.log(10), tot[i]),
                xytext=(25, 10), textcoords="offset points", fontsize=7.5)
    a1.set_yscale("symlog", linthresh=1e4)
    a1.set_xlabel(r"$\log_{10}(1/R\,/\,{\rm GeV})$")
    a1.set_ylabel("action")
    a1.set_title(r"(a) smooth model $\lambda=\lambda_0+\frac{1}{2}\lambda_2(x-x_0)^2$, $q=25$")
    a1.legend(fontsize=7.5, loc="upper right")
    a1.grid(alpha=0.3, lw=0.3)
    if legacy is not None:
        d = legacy[(legacy.legacy_S_exact != -1) & legacy.S_num.notna()].copy()
        neg = d.legacy_S_exact < 0
        sent = d.legacy_S_exact >= 1e99
        ok = ~neg & ~sent
        agree = d.legacy_status == d.status
        a2.semilogx(d.S_num[ok & agree], (d.legacy_S_exact / d.S_num - 1)[ok & agree], "o", ms=3, color="0.4",
                    label="same classification")
        a2.semilogx(d.S_num[ok & ~agree], (d.legacy_S_exact / d.S_num - 1)[ok & ~agree], "o", ms=4, color="#b2182b",
                    label="different classification")
        if neg.any():
            a2.semilogx(d.S_num[neg], np.full(neg.sum(), -0.95), "v", ms=6, color="k", label=f"legacy S < 0 ({int(neg.sum())})")
        a2.axhline(0, color="k", lw=0.6)
        a2.set_xlabel(r"$S_{\rm num}$ (new reference, cw1)")
        a2.set_ylabel(r"$S_{\rm exact}^{\rm legacy}/S_{\rm num}-1$")
        a2.set_title("(b) legacy vs new on the controlled region (0.5 GeV grid)")
        a2.set_ylim(-1.0, 1.2)
        a2.legend(fontsize=7.5)
        a2.grid(alpha=0.3, lw=0.3, which="both")
    fig.savefig(os.path.join(FIG, "legacy_failure.png"))
    plt.close(fig)


def fig_near_critical(model, smlines, scans):
    """Approach to the absolute-stability boundary: (a) S/S_FL vs q, (b) errors of the estimates vs q,
    (c) actions along fixed-M_h lines vs delta = M_t - M_t^crit."""
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    a1, a2, a3 = axs
    qq = np.logspace(-1.4, 3.2, 300)
    # (a) universality of S/S_FL as a function of q
    for name, df in scans.items():
        d = df[df.S_num.notna() & (df.lambda2 > 0) & (df.planck_sensitive == 0)]
        a1.loglog(d.lambda_min.abs() / d.lambda2, d.S_num / d.S_FL - 1, ".", ms=1.5, alpha=0.25,
                  color="0.55" if name == "cw1" else "0.75", rasterized=True)
    m = model.sort_values("q")
    a1.loglog(m.q, m.S_num / m.S_FL - 1, "k-o", ms=4, lw=1.4, label="parabolic model (numerical)")
    a1.loglog(qq, (13 / 72) / qq, "k--", lw=0.9, label=r"thick-wall law $\frac{13}{72}q^{-1}$")
    qt = np.logspace(math.log10(0.02), math.log10(0.2), 60)       # leading-order Coleman thin-wall action
    tw = (27 * 64 / (2 * 18 ** 4)) * 3 / 8 / qt ** 2 - 1            # S_TW/S_FL - 1 for the parabolic model
    a1.loglog(qt[tw > 0], tw[tw > 0], "k:", lw=1.1, label=r"leading thin-wall $S_{\rm TW}/S_{\rm FL}-1$")
    styles = {("cw1", 120.0): ("#1b9e77", "s"), ("cw1", 125.15): ("#d95f02", "o"), ("cw1", 128.0): ("#7570b3", "^"),
              ("tree", 125.15): ("#e7298a", "D")}
    for (presc, Mh), sub in smlines.groupby(["prescription", "Mh"]):
        c, mk = styles.get((presc, Mh), ("k", "x"))
        s = sub[sub.S_num.notna() & (sub.lambda2 > 0)]
        q = s.lambda_min.abs() / s.lambda2
        ps = s.planck_sensitive.astype(bool)
        a1.loglog(q[~ps], (s.S_num / s.S_FL - 1)[~ps], mk, color=c, ms=5, label=f"SM {presc}, $M_h$={Mh:g}")
        a1.loglog(q[ps], (s.S_num / s.S_FL - 1)[ps], mk, mfc="none", color=c, ms=5)
    a1.set_xlabel(r"$q=|\lambda_{\min}|/\lambda''_{\min}$")
    a1.set_ylabel(r"$S_{\rm num}/S_{\rm FL}-1$")
    a1.set_title("(a) excess of the true action over the FL bound")
    a1.legend(fontsize=7, loc="upper right")
    a1.grid(alpha=0.3, which="both", lw=0.3)
    # (b) errors of the estimates
    m = m.assign(S_E1=m.S_FL * (1 + (13 / 72) / m.q))                       # Espinosa (2020) eq. (3.21)
    smlines = smlines.assign(S_E1=smlines.S_FL * (1 + (13 / 72) * smlines.lambda2 / smlines.lambda_min.abs()))
    for col, c, lab in (("S_FL", "#b2182b", r"$S_{\rm FL}$ (lower bound)"),
                        ("S_E1", "#e6ab02", r"$S_{\rm E1}$ (first order, Espinosa 2020)"),
                        ("S_curv", "#2166ac", r"$S_{\rm curv}$ (resummed)"),
                        ("S_CGM", "#1b9e77", r"$S_{\rm CGM}$ (FL trial)"), ("S_CGM2", "#7570b3", r"$S_{\rm CGM2}$ (two-parameter trial)")):
        err = (m[col] - m.S_num).abs() / m.S_num
        a2.loglog(m.q, err.where(np.isfinite(m[col])), "-o", color=c, ms=3, lw=1.3, label=lab + " (model)")
        s = smlines[smlines.S_num.notna() & (smlines.lambda2 > 0)]
        q = s.lambda_min.abs() / s.lambda2
        a2.loglog(q, ((s[col] - s.S_num).abs() / s.S_num).where(np.isfinite(s[col])), ".", color=c, ms=4, alpha=0.7)
    a2.axvline(13 / 72, color="#1b9e77", ls=":", lw=1)
    a2.text(13 / 72 * 1.06, 2e-9, r"$q=13/72$: $S_{\rm CGM}=\infty$ below", rotation=90, fontsize=7, color="#1b9e77")
    a2.axvline(3 / 32, color="#7570b3", ls=":", lw=1)
    a2.text(3 / 32 * 1.06, 2e-9, r"$q=3/32$: $S_{\rm CGM2}=\infty$ below", rotation=90, fontsize=7, color="#7570b3")
    a2.set_xlabel(r"$q=|\lambda_{\min}|/\lambda''_{\min}$")
    a2.set_ylabel("|estimate - numerical| / numerical")
    a2.set_title("(b) accuracy of the estimates (lines: model; dots: SM lines)")
    a2.set_ylim(1e-10, 300)
    a2.plot([], [], " ", label=r"(model: $S_{\rm curv}=S_{\rm CGM}$ exactly; lines coincide)")
    a2.legend(fontsize=7, loc="upper right")
    a2.grid(alpha=0.3, which="both", lw=0.3)
    # (c) along the SM lines
    for (presc, Mh), sub in smlines.groupby(["prescription", "Mh"]):
        c, mk = styles.get((presc, Mh), ("k", "x"))
        s = sub.sort_values("delta")
        good = s.S_num.notna()
        a3.loglog(s.delta[good], s.S_num[good], mk + "-", color=c, ms=4, lw=1, label=f"{presc}, $M_h$={Mh:g}: numerical")
        a3.loglog(s.delta, s.S_FL, "--", color=c, lw=0.8)
        ps = s.planck_sensitive.astype(bool) & good
        a3.loglog(s.delta[ps], s.S_num[ps], "o", mfc="none", mec="k", ms=8)
        nfb = s.status == "no-finite-bounce"
        a3.loglog(s.delta[nfb], s.S_FL[nfb], "x", color=c, ms=6, mew=1.2)
        fail = ~good & ~nfb
        a3.loglog(s.delta[fail], s.S_FL[fail], "+", color=c, ms=8, mew=1.2)
    a3.plot([], [], "k--", lw=0.8, label=r"$S_{\rm FL}$ (rigorous lower bound)")
    a3.plot([], [], "o", mfc="none", mec="k", label="Planck-sensitive")
    a3.plot([], [], "kx", label=r"no finite bounce (Planck plateau), at $S_{\rm FL}$")
    a3.plot([], [], "k+", ms=8, label=r"unresolved ($q\lesssim0.05$), at $S_{\rm FL}$")
    a3.set_xlabel(r"$\delta=M_t-M_t^{\rm crit}$ [GeV]")
    a3.set_ylabel("bounce action S")
    a3.set_title("(c) approach to the stability boundary")
    a3.legend(fontsize=7, loc="upper right")
    a3.grid(alpha=0.3, which="both", lw=0.3)
    fig.savefig(os.path.join(FIG, "near_critical.png"))
    plt.close(fig)


def fig_variational_failure(ncs):
    """Mechanism of the variational failure at a near-critical SM point: (a) numerical bounce vs the optimal
    one-parameter (FL) and two-parameter trial profiles; (b) distribution of ln phi under the weight phi^4 r^3 dr,
    which decides the sign of V (the estimate is finite only if this weight sits mostly where lambda < 0)."""
    sys.path.insert(0, os.path.join(ROOT, "python"))
    from smvac.reference import analytic as A, bounce as B, potential as P, sm
    row = ncs[(ncs.prescription == "cw1") & np.isclose(ncs.Mh, 125.15) & ncs.S_num.notna()]
    row = row.iloc[(row.lambda_min.abs() / row.lambda2 - 0.35).abs().argmin()]
    pot = P.sm_potential(sm.Inputs(125.15, float(row.Mt)), "cw1")
    an = A.estimates(pot)
    q = abs(an.lambda_min) / an.lambda2
    r = None
    for N in (256, 512, 1024):
        r = B.solve_m4(pot, an.lambda_min, an.lnh_CGM, N=N)
        if r.ok and abs(r.derrick) < 1e-9 and abs(r.dilatation) < 1e-9:
            break
    from scipy.interpolate import BarycentricInterpolator
    tn, fn = r.profile
    tau = np.linspace(tn.min(), tn.max(), 40001)
    curves = {"numerical bounce (M4)": (tau, BarycentricInterpolator(tn, fn)(tau), r.S, "#000000", "-")}
    S2, lnh2, n2 = A.cgm2_estimate(pot, an.lnh_CGM)
    for lab, lnh, n, S, c in ((r"optimal FL trial ($S_{\rm CGM}$)", an.lnh_CGM, 1.0, an.S_CGM, "#1b9e77"),
                              (rf"optimal two-parameter trial, $n={n2:.2f}$ ($S_{{\rm CGM2}}$)", lnh2, n2, S2, "#7570b3")):
        T, V = A.cgm2_trial(pot, lnh, n)
        lnR = 0.5 * math.log(T / (2.0 * abs(V)))
        t = np.linspace(lnR - 15, lnR + 15, 40001)
        lw = np.array([A._log_w(u, n) for u in t - lnR])
        curves[lab] = (t, np.exp(t + lnh + lw), S, c, "--")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    x_axis = np.linspace(an.x_min - 4, an.x_min + 1.5, 400)
    bins = np.linspace(x_axis[0], x_axis[-1], 111)
    for lab, (t, f, S, c, ls) in curves.items():
        rr, phi = np.exp(t), f * np.exp(-t)
        err = "" if lab.startswith("numerical") else f", {100 * (S / r.S - 1):+.1f}%"
        a1.loglog(rr, phi, ls, color=c, lw=1.6, label=f"{lab}: S = {S:.5g}{err}")
        x = np.log(np.abs(f)) - t
        wgt = f ** 4 * (t[1] - t[0])
        hist, _ = np.histogram(x, bins=bins, weights=wgt)
        hist = hist / hist.sum() / (bins[1] - bins[0])
        mean = np.sum(wgt * x) / np.sum(wgt)
        hv = 0.5 * np.sum(wgt * (x - mean) ** 2) / np.sum(wgt)
        a2.step(0.5 * (bins[1:] + bins[:-1]) / math.log(10), hist, where="mid", color=c, ls=ls, lw=1.4,
                label=rf"{lab.split(' (')[0]}: $\frac{{1}}{{2}}{{\rm Var}}(\ln\phi)={hv:.3f}$")
    a1.axhline(math.exp(an.x_min), color="#b2182b", lw=0.8, ls=":", label=r"field at $\lambda_{\min}$")
    a1.set_xlabel(r"$r$ [GeV$^{-1}$]")
    a1.set_ylabel(r"$\phi(r)$ [GeV]")
    phi0 = math.exp(r.lnphi0)
    a1.set_ylim(phi0 * 1e-4, phi0 * 3)
    a1.set_xlim(math.exp(r.lnR) * 10 ** -2.5, math.exp(r.lnR) * 10 ** 2.0)
    a1.set_title(rf"(a) profiles at $M_h=125.15$, $M_t=M_t^{{\rm crit}}+{row.delta:.4f}$ GeV ($q={q:.3f}$)")
    a1.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    a1.grid(alpha=0.3, which="both", lw=0.3)
    a2.set_xlabel(r"$\log_{10}(\phi/{\rm GeV})$")
    a2.set_ylabel(r"weight $\phi^4r^3\,dr$ per unit $\ln\phi$ (normalised)")
    a2.set_title(r"(b) where each profile puts its potential-energy weight")
    ax2 = a2.twinx()
    lam = np.array([pot.ev(x)[0] for x in x_axis])
    ax2.plot(x_axis / math.log(10), lam, color="#b2182b", lw=1.0)
    ax2.axhline(0, color="#b2182b", lw=0.5, ls=":")
    ax2.fill_between(x_axis / math.log(10), lam, 0, where=lam < 0, color="#b2182b", alpha=0.15)
    ax2.set_ylabel(r"$\lambda_{\rm eff}(\phi)$ (red; shaded: $\lambda<0$)", color="#b2182b")
    ax2.set_ylim(-4 * abs(an.lambda_min), 12 * abs(an.lambda_min))
    a2.legend(fontsize=7, loc="upper left")
    a2.text(0.02, 0.30, rf"$V<0$ requires the weight to sit mostly where $\lambda<0$;" "\n"
            rf"for quadratic $\lambda$: $\frac{{1}}{{2}}\langle(\ln\phi-x_{{\min}})^2\rangle<q={q:.3f}$",
            transform=a2.transAxes, ha="left", va="bottom", fontsize=7, bbox=dict(fc="w", ec="0.7", lw=0.5))
    a2.grid(alpha=0.3, lw=0.3)
    fig.savefig(os.path.join(FIG, "variational_failure.png"))
    plt.close(fig)


def main():
    cw1, meta = read("results/controlled_cw1.csv")
    tree, _ = read("results/controlled_tree.csv")
    comp, _ = read("validation/results/method_comparison.csv")
    conv, _ = read("validation/results/convergence.csv")
    expl, _ = read("results/exploratory_cw1.csv")
    legacy, _ = read("validation/results/legacy_comparison.csv")
    fb, _ = read("results/controlled_cw1_fallback.csv")
    if cw1 is not None:
        cw1 = apply_fallback(cw1, fb)
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
    fig_legacy(legacy)
    ncm, _ = read("validation/results/near_critical_model.csv")
    ncs, _ = read("validation/results/near_critical_sm.csv")
    if ncm is not None and ncs is not None and cw1 is not None:
        fig_near_critical(ncm, ncs, {k: v for k, v in (("cw1", cw1), ("tree", tree)) if v is not None})
    if ncs is not None:
        fig_variational_failure(ncs)
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
