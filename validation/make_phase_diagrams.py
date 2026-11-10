"""Phase diagrams in the style of the original SMVac figures, computed with the validated numerical bounce.

Every figure shows the measured (M_h, M_t) with its 1 sigma and 2 sigma regions, drawn as the joint confidence
ellipses of two independent Gaussian measurements (grey areas, as in Buttazzo et al. 2013, fig. 3).

Figures (figures/):
  phase_diagram_full.png                    [0,250]^2 GeV: stability regions and boundaries; inset: magnified view
                                            around the measured point
  phase_diagram_full_contours.png           the same with contours of constant bounce action S
  phase_diagram_full_instability_scale.png  the same with contours of the instability scale Lambda_I
                                            (as in Buttazzo et al. 2013, fig. 3 left)
  phase_diagram_closeup.png                 [110,140]x[160,185] GeV, the window of the original close-ups
  phase_diagram_closeup_contours.png        the close-up with contours of constant S
  phase_diagram_zoom.png                    controlled region [120,130]x[166,182] GeV with the instability scale
                                            (as in Buttazzo et al. 2013, fig. 3 right); from the close-up scan if
                                            present, otherwise from the controlled scan
  overlay_full.png, overlay_closeup.png     analytical (Fubini-Lipatov, literature scale) vs numerical classification;
                                            disagreements in blue
  approximation_errors_full.png,            where the closed-form actions fail: relative errors of S_FL, S_E1, S_curv
  approximation_errors_closeup.png          and S_CGM along M_h ~ 125.15 GeV versus the distance from the stability
                                            boundary, next to a map of the S_FL error (full plane / close-up)
Inputs (CSV, optionally .gz): results/full_plane_cw1.csv, results/closeup_cw1.csv, results/controlled_cw1.csv and,
where present, the shooting-fallback results <scan>_fallback.csv (validation/fallback_shooting.py); the near-critical
line validation/results/near_critical_sm.csv extends the cut of the error figures to M_t - M_t^crit = 1e-4 GeV.

Points that neither solver resolves are classified by the action bound where it decides alone (hatched; see
apply_action_bound): S >= 8 pi^2/(3|lambda_min|) and the decay probability needs S < S_CRIT_PLANCK even for a bounce
at the Planck scale. Points the bound cannot decide remain 'unresolved'.
Numbers (results/phase_diagrams_summary.json, quoted in README section 10.7): status counts, fallback and action-bound
counts, agreement of the analytical classification, S_FL errors, the measured point and its sigma regions, and the
boundaries at M_h = 125.2 GeV.
Usage:  python validation/make_phase_diagrams.py [--full CSV] [--closeup CSV] [--zoom CSV] [--outdir DIR]
        [--summary JSON]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "validation"))
import fallback_shooting as FS  # noqa: E402  (H0 and M_Pl conventions of the decay probability)
import make_figures as MF  # noqa: E402  (read, grid, apply_fallback; sets the Agg backend)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import PatchCollection  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm, to_rgba  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Ellipse, Patch, Rectangle  # noqa: E402

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11, "axes.labelsize": 13,
                     "axes.titlesize": 13, "legend.fontsize": 9, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True, "savefig.dpi": 200, "savefig.bbox": "tight",
                     "hatch.linewidth": 0.7})

# Measured masses: Review of Particle Physics, 2026 edition, summary tables (M_t from direct measurements).
MEASURED = {"Mh": 125.13, "dMh": 0.11, "Mt": 172.60, "dMt": 0.27, "source": "PDG 2026"}
# n-sigma joint region of two independent Gaussian measurements: chi^2 <= -2 ln(1 - erf(n/sqrt 2)), i.e.
# Delta chi^2 = 2.30 (68.3%) and 6.18 (95.4%).
SIGMA = {n: -2.0 * math.log(1.0 - math.erf(n / math.sqrt(2.0))) for n in (1, 2)}

STATUS = ["stable", "metastable", "metastable-bound", "unstable", "non-perturbative", "no-bounce", "solver-failure",
          "invalid-input"]
COLORS = {"stable": "#78f078", "metastable": "#ffff66", "metastable-bound": "#ffff66", "unstable": "#ff7878",
          "non-perturbative": "#3c3c3c", "no-bounce": "#c8c8c8", "solver-failure": "#b000b0", "invalid-input": "#ffffff",
          "disagreement": "#1f3cff"}
HATCHES = {"metastable-bound": "////"}    # drawn over the category colour
HATCH_COLOR = "0.3"
LABELS = {"stable": "Stability", "metastable": "Meta-stability",
          "metastable-bound": r"Meta-stability from the action bound $S\geq8\pi^2/(3|\lambda_{\min}|)$ (bounce not "
                              r"resolved)",
          "unstable": "Instability", "non-perturbative": "Non-perturbative",
          "no-bounce": r"Planck-scale dominated (no finite bounce, $\lambda_{\min}$ at $M_{\rm Pl}$)",
          "solver-failure": "Unresolved (solver failure)", "invalid-input": r"Invalid input ($y_t(M_t)\leq0$)"}
BOUNCE = ["metastable", "unstable"]
LAMBDA_VALID = ["stable"] + BOUNCE + ["metastable-bound", "no-bounce", "solver-failure"]   # lambda(mu) computed
# Action at which p = 1 for a bounce at the Planck scale (conventions of validation/fallback_shooting.py): a larger
# action gives p < 1 for every bounce scale up to M_Pl.
S_CRIT_PLANCK = 4.0 * (FS.LN_M_PLANCK - math.log(FS.H0_GEV)) + math.log(0.15)
XLABEL, YLABEL = r"Higgs pole mass $M_h$ in GeV", r"Top pole mass $M_t$ in GeV"
CONTROLLED = (120.0, 130.0, 166.0, 182.0)


def load(path):
    """Scan output with the accepted shooting-fallback results applied (if <stem>_fallback.csv exists)."""
    df, meta = MF.read(path)
    if df is None:
        return None, {}
    stem = path[:-3] if path.endswith(".gz") else path
    stem = stem[:-4] if stem.endswith(".csv") else stem
    fb, _ = MF.read(stem + "_fallback.csv")
    df = apply_action_bound(MF.apply_fallback(df, fb), fb)
    meta["fallback_points"] = int(df.fallback.sum())
    meta["bound_points"] = int(df.bound.sum())
    meta["step"] = float(np.round(np.diff(np.sort(df.Mh.unique())).min(), 6))
    return df, meta


def apply_action_bound(df, fb):
    """Solver failures that are metastable by the action bound alone.

    V(h) >= lambda_min h^4/4 for every h and the bounce action grows with the potential, so S >= S_FL =
    8 pi^2/(3|lambda_min|) (Fubini-Lipatov at lambda_min; every computed bounce in our scans satisfies it). Where
    S_FL > S_CRIT_PLANCK, p < 1 for every bounce scale up to M_Pl. Applied only to an interior lambda minimum (the
    fallback found the Coleman-Glaser-Martin optimum at finite field), so the Planck-plateau case keeps its category."""
    df = df.copy()
    df["bound"] = 0
    if fb is None or "CGM_interior" not in fb:
        return df
    interior = {(round(r.Mh, 6), round(r.Mt, 6)) for r in fb.itertuples() if r.CGM_interior == 1}
    cand = df.index[(df.status == "solver-failure") & (df.lambda_min < 0) & (df.S_FL > S_CRIT_PLANCK)]
    hit = [j for j in cand if (round(df.at[j, "Mh"], 6), round(df.at[j, "Mt"], 6)) in interior]
    df.loc[hit, ["status", "bound"]] = ["metastable-bound", 1]
    return df


def masked(df, col, keep):
    """Grid of `col`, masked outside the statuses in `keep` (fields such as lambda_min are set to 0 there)."""
    mh, mt, z = MF.grid(df.assign(_v=df[col].where(df.status.isin(keep))), "_v")
    return mh, mt, np.ma.masked_invalid(z)


def draw_categories(ax, df, column, order):
    codes = {s: i for i, s in enumerate(order)}
    mh, mt, z = MF.grid(df.assign(_c=df[column].map(codes)), "_c")
    cmap = ListedColormap([COLORS[s] for s in order])
    ax.pcolormesh(mh, mt, z, cmap=cmap, norm=BoundaryNorm(np.arange(-0.5, len(order) + 0.5), cmap.N),
                  shading="nearest", rasterized=True)
    dx, dy = np.diff(mh).min(), np.diff(mt).min()
    for s, hatch in HATCHES.items():        # hatching over the cells of these categories (aligned across cells)
        sel = df[df[column] == s]
        if s in codes and len(sel):
            cells = [Rectangle((x - dx / 2, y - dy / 2), dx, dy) for x, y in zip(sel.Mh, sel.Mt)]
            ax.add_collection(PatchCollection(cells, facecolor="none", edgecolor=HATCH_COLOR, linewidth=0,
                                              hatch=hatch))


def draw_boundaries(ax, df, lw=1.1):
    """lambda_min = 0 (absolute stability) and p = 1 (numerical tunnelling probability over the past light cone)."""
    mh, mt, lam = masked(df, "lambda_min", LAMBDA_VALID)
    ax.contour(mh, mt, lam, levels=[0.0], colors="k", linewidths=lw)
    mh, mt, p = masked(df, "log10_p_num", BOUNCE)
    ax.contour(mh, mt, p, levels=[0.0], colors="k", linewidths=lw)


def draw_measured(ax, star_size=13, lw=1.0):
    m = MEASURED
    for n, fc in ((2, to_rgba("0.45", 0.40)), (1, to_rgba("0.20", 0.60))):
        a = math.sqrt(SIGMA[n])
        ax.add_patch(Ellipse((m["Mh"], m["Mt"]), 2 * a * m["dMh"], 2 * a * m["dMt"], fc=fc, ec="k", lw=lw, zorder=5))
    ax.plot(m["Mh"], m["Mt"], marker="*", ms=star_size, mfc="#1f3cff", mec="k", mew=0.8, ls="none", zorder=6)


def measured_handles(star_size=13):
    m = MEASURED
    return ([Line2D([], [], marker="*", ms=star_size, mfc="#1f3cff", mec="k", mew=0.8, ls="none"),
             Patch(fc=to_rgba("0.20", 0.60), ec="k", lw=1.0), Patch(fc=to_rgba("0.45", 0.40), ec="k", lw=1.0)],
            [rf"{m['source']}: $M_h={m['Mh']:.2f}\pm{m['dMh']:.2f}$, $M_t={m['Mt']:.2f}\pm{m['dMt']:.2f}$ GeV",
             rf"$1\sigma$ region (68.3% CL, $\Delta\chi^2={SIGMA[1]:.2f}$)",
             rf"$2\sigma$ region (95.4% CL, $\Delta\chi^2={SIGMA[2]:.2f}$)"])


def region_handles(df, column="status", order=STATUS, labels=LABELS):
    present = [s for s in order if (df[column] == s).any()]
    return ([Patch(fc=COLORS[s], ec="k", lw=0.4, hatch=HATCHES.get(s)) for s in present],
            [labels[s] for s in present])


def finish(ax, xlim, ylim, title, subtitle):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_title(title, pad=8 + 11 * (subtitle.count("\n") + 1))
    ax.text(0.0, 1.008, subtitle, transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5, color="0.25",
            linespacing=1.3)
    ax.grid(True, linestyle="--", alpha=0.45, lw=0.6)


def subtitle(meta, extra="", split=False):
    """Provenance line(s) above a panel; split=True puts each clause on its own line (narrow panels)."""
    parts = [f"numerical bounce (Chebyshev collocation) on a {meta['step']:g} GeV grid",
             f"prescription {meta.get('prescription', '?')}, matching {meta.get('matching', '?')}"]
    counts = []
    if meta.get("fallback_points"):
        counts.append(f"{meta['fallback_points']} points from the shooting fallback")
    if meta.get("bound_points"):
        counts.append(f"{meta['bound_points']} {'' if counts else 'points '}from the action bound")
    if counts:
        parts.append(", ".join(counts))
    return ("\n" if split else "; ").join(parts) + ("\n" + extra if extra else "")


INSET = [0.73, 0.44, 0.25, 0.26]     # axes fraction, over the non-perturbative region, clear of the legend


def add_inset(ax, fine, window, bounds=INSET, column="status", order=STATUS, boundaries=True):
    """Magnified view of `window` (Mh0, Mh1, Mt0, Mt1) drawn from `fine`, placed at `bounds` (axes fraction) on a
    white card that also holds the tick labels."""
    x, y, w, h = bounds
    ax.add_patch(plt.Rectangle((x - 0.052, y - 0.058), w + 0.062, h + 0.068, transform=ax.transAxes, fc="w", ec="k",
                               lw=0.6, zorder=4))
    ins = ax.inset_axes(bounds)
    draw_categories(ins, fine, column, order)
    if boundaries:
        draw_boundaries(ins, fine, lw=0.9)
    draw_measured(ins, star_size=9)
    ins.set_xlim(window[0], window[1])
    ins.set_ylim(window[2], window[3])
    ins.tick_params(labelsize=8)
    ins.set_facecolor("w")
    ins.grid(True, linestyle="--", alpha=0.45, lw=0.5)
    ins.text(0.02, 0.98, f"magnified ({fine.attrs.get('step', '?'):g} GeV grid)", transform=ins.transAxes, ha="left",
             va="top", fontsize=7.5, bbox=dict(fc="w", ec="none", alpha=0.8, pad=1))
    ax.indicate_inset_zoom(ins, edgecolor="k", lw=0.8, alpha=0.9)
    return ins


def labelled_contours(ax, draw, **kw):
    """Contours from draw() with inline labels at matplotlib's automatic positions, except that a label overlapping
    the measured point or its sigma regions moves along its contour to the nearest position clear of them."""
    cs = draw()
    trial = ax.clabel(cs, inline=False, **kw)
    renderer = ax.figure.canvas.get_renderer()
    to_data = ax.transData.inverted()
    a = math.sqrt(SIGMA[2])
    ex, ey = a * MEASURED["dMh"], a * MEASURED["dMt"]           # half-widths of the 2 sigma region
    (xlo, xhi), (ylo, yhi) = ax.get_xlim(), ax.get_ylim()
    paths = [p.vertices for p in cs.get_paths()]
    positions = []
    for t in trial:
        x, y = t.get_position()
        (x0, y0), (x1, y1) = to_data.transform(t.get_window_extent(renderer).get_points())
        hw, hh = abs(x1 - x0) / 2, abs(y1 - y0) / 2            # half-size of the label
        clear = lambda px, py: ((np.abs(px - MEASURED["Mh"]) > ex + 3 * hw)   # a label width of space beside
                                | (np.abs(py - MEASURED["Mt"]) > ey + 3 * hh))  # or a label height above/below
        if not clear(x, y):
            dist = lambda v: np.hypot((v[:, 0] - x) / (xhi - xlo), (v[:, 1] - y) / (yhi - ylo))
            v = min((v for v in paths if len(v)), key=lambda v: dist(v).min())   # this label's contour
            ok = v[clear(v[:, 0], v[:, 1]) & (v[:, 0] > xlo + hw) & (v[:, 0] < xhi - hw)
                   & (v[:, 1] > ylo + hh) & (v[:, 1] < yhi - hh)]
            if len(ok):
                x, y = ok[np.argmin(dist(ok))]
        positions.append((x, y))
    cs.remove()                                                 # with the trial labels
    cs = draw()
    ax.clabel(cs, manual=positions, inline=True, **kw)
    return cs


def action_contours(ax, df, levels, color="#0b2a9e", fontsize=8):
    """Contours of constant S; drawn on 1/S (proportional to |lambda_min| at leading order, so it interpolates
    smoothly between grid points)."""
    mh, mt, s = masked(df, "S_num", BOUNCE)
    inv = np.ma.masked_where(~(s > 0), 1.0 / s)
    return labelled_contours(
        ax, lambda: ax.contour(mh, mt, inv, levels=sorted(1.0 / np.asarray(levels)), colors=color, linewidths=0.9),
        fmt=lambda v: f"S={1.0 / v:,.0f}", fontsize=fontsize, inline_spacing=3)


def instability_contours(ax, df, levels):
    mh, mt, li = masked(df, "log10_mu_inst", [s for s in LAMBDA_VALID if s != "stable"])
    return labelled_contours(
        ax, lambda: ax.contour(mh, mt, li, levels=levels, colors="k", linewidths=0.8, linestyles=":"),
        fmt=lambda v: rf"$10^{{{v:.0f}}}$", fontsize=8)


def legend(ax, handles, labels, loc, ncol=1, fontsize=9):
    if loc == "below":
        leg = ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, fontsize=fontsize,
                        framealpha=0.95, edgecolor="k")
        leg.set_zorder(20)
        return leg
    leg = ax.legend(handles, labels, loc=loc, ncol=ncol, fontsize=fontsize, framealpha=0.95, edgecolor="k")
    leg.set_zorder(20)
    return leg


BOUNDARY_HANDLE = (Line2D([], [], color="k", lw=1.1), r"Boundaries: $\lambda_{\min}=0$ and $\wp=1$")
EXTRAPOLATION = "matching extrapolated away from the measured masses (controlled region 120-130 x 166-182 GeV)"


def inset_window():
    m = MEASURED
    return (m["Mh"] - 2.6, m["Mh"] + 2.6, m["Mt"] - 3.4, m["Mt"] + 3.0)


def fig_full(df, meta, fine, out, contours=None):
    fig, ax = plt.subplots(figsize=(10, 8.2))
    draw_categories(ax, df, "status", STATUS)
    draw_boundaries(ax, df)
    draw_measured(ax)
    handles, labels = region_handles(df)
    title = "SM vacuum phase diagram"
    if contours == "action":
        action_contours(ax, df, [75, 100, 150, 200, 300, 450, 1000, 3000, 10000, 100000])
        handles.append(Line2D([], [], color="#0b2a9e", lw=0.9))
        labels.append("Contours of constant bounce action $S$")
        title += " with contours of constant action"
    elif contours == "instability":
        instability_contours(ax, df, [3, 4, 5, 6, 8, 10, 12, 14, 16, 18])
        handles.append(Line2D([], [], color="k", lw=0.8, ls=":"))
        labels.append(r"Instability scale $\Lambda_I$ in GeV")
        title += r" with the instability scale $\Lambda_I$"
    handles.append(BOUNDARY_HANDLE[0])
    labels.append(BOUNDARY_HANDLE[1])
    h2, l2 = measured_handles()
    legend(ax, handles + h2, labels + l2, "lower right", fontsize=8)
    add_inset(ax, fine, inset_window())
    finish(ax, (0, 250), (0, 250), title, subtitle(meta, EXTRAPOLATION))
    fig.savefig(out)
    plt.close(fig)


def fig_window(df, meta, window, out, title, contours=None, levels=None, star_size=11, loc="lower right"):
    # only the window, so that contour labels are placed where they are visible and the counts refer to it
    df = df[(df.Mh >= window[0]) & (df.Mh <= window[1]) & (df.Mt >= window[2]) & (df.Mt <= window[3])]
    meta = dict(meta, fallback_points=int(df.fallback.sum()), bound_points=int(df.bound.sum()))
    fig, ax = plt.subplots(figsize=(10, 8.2))
    draw_categories(ax, df, "status", STATUS)
    draw_boundaries(ax, df, lw=1.3)
    handles, labels = region_handles(df)
    if contours == "action":
        action_contours(ax, df, levels, fontsize=8.5)
        handles.append(Line2D([], [], color="#0b2a9e", lw=0.9))
        labels.append("Contours of constant bounce action $S$")
    elif contours == "instability":
        instability_contours(ax, df, levels)
        handles.append(Line2D([], [], color="k", lw=0.8, ls=":"))
        labels.append(r"Instability scale $\Lambda_I$ in GeV")
    draw_measured(ax, star_size=star_size)
    handles.append(BOUNDARY_HANDLE[0])
    labels.append(BOUNDARY_HANDLE[1])
    h2, l2 = measured_handles(star_size)
    legend(ax, handles + h2, labels + l2, loc)
    finish(ax, window[:2], window[2:], title, subtitle(meta))
    fig.savefig(out)
    plt.close(fig)


def overlay_status(df):
    """Analytical classification: S_FL with Lambda_B = argmin lambda (the literature convention); it differs from
    the numerical one only for points with a bounce (stability and perturbativity depend on lambda(mu) alone)."""
    num = df.status
    ana = num.copy()
    b = num.isin(BOUNCE)
    ana[b] = np.where(df.loc[b, "log10_p_FL_lit"] < 0, "metastable", "unstable")
    return num.where(~b | (num == ana), "disagreement")


def fig_overlay(df, meta, window, out, title, fine_window=False, fine=None):
    d = df.assign(overlay=overlay_status(df))
    order = STATUS + ["disagreement"]
    labels = dict(LABELS)
    labels.update({"stable": "Stability (both agree)", "metastable": "Meta-stability (both agree)",
                   "unstable": "Instability (both agree)",
                   "disagreement": r"Disagreement: $S_{\rm FL}$ (literature scale) vs numerical bounce"})
    fig, ax = plt.subplots(figsize=(10, 8.2))
    draw_categories(ax, d, "overlay", order)
    draw_measured(ax, star_size=11 if fine_window else 13)
    sel = d[(d.Mh >= window[0]) & (d.Mh <= window[1]) & (d.Mt >= window[2]) & (d.Mt <= window[3])]
    handles, labs = region_handles(sel, "overlay", order, labels)
    n_dis = int((sel.overlay == "disagreement").sum())
    n_b = int(sel.status.isin(BOUNCE).sum())
    h2, l2 = measured_handles(11 if fine_window else 13)
    legend(ax, handles + h2, labs + l2, "lower right", fontsize=9 if fine_window else 8)
    note = f"{n_dis} of {n_b} bounce points classified differently"
    if fine is not None:
        f = fine.assign(overlay=overlay_status(fine))
        f.attrs["step"] = fine.attrs.get("step")
        add_inset(ax, f, inset_window(), column="overlay", order=order, boundaries=False)
        note += "; " + EXTRAPOLATION
    finish(ax, window[:2], window[2:], title, subtitle(meta, note))
    fig.savefig(out)
    plt.close(fig)


ESTIMATES = {"S_FL": ("#b2182b", r"$S_{\rm FL}=8\pi^2/(3|\lambda_{\min}|)$ (lower bound)"),   # colours as in
             "S_E1": ("#e6ab02", r"$S_{\rm E1}=S_{\rm FL}(1+\frac{13}{72q})$ (first order)"),  # make_figures.py
             "S_curv": ("#2166ac", r"$S_{\rm curv}$ (resummed)"),
             "S_CGM": ("#1b9e77", r"$S_{\rm CGM}$ (variational)")}
NEAR_CRITICAL = os.path.join(ROOT, "validation", "results", "near_critical_sm.csv")


def with_e1(d):
    """Adds S_E1 = S_FL (1 + 13/(72 q)), q = |lambda_min| / lambda''_min (Espinosa 2020, eq. 3.21)."""
    return d.assign(S_E1=d.S_FL * (1 + 13 / (72 * d.lambda_min.abs() / d.lambda2)))


def critical_cut(scan):
    """Bounce points of the M_h line of `scan` closest to the measured M_h, with delta = M_t - M_t^crit (lambda_min = 0
    interpolated linearly), and (delta, S) where p = 1 on that line (None if the line does not reach it)."""
    mhs = np.sort(scan.Mh.unique())
    mh0 = mhs[np.argmin(np.abs(mhs - MEASURED["Mh"]))]
    line = scan[np.isclose(scan.Mh, mh0) & scan.status.isin(LAMBDA_VALID)].sort_values("Mt")
    lam, mt = line.lambda_min.to_numpy(), line.Mt.to_numpy()
    k = np.nonzero((lam[:-1] > 0) & (lam[1:] <= 0))[0][0]
    mt_crit = mt[k] + lam[k] / (lam[k] - lam[k + 1]) * (mt[k + 1] - mt[k])
    b = with_e1(line[line.status.isin(BOUNCE) & (line.S_num > 0)]).assign(delta=lambda d: d.Mt - mt_crit)
    lp, dl = b.log10_p_num.to_numpy(), b.delta.to_numpy()
    j = np.nonzero((lp[:-1] < 0) & (lp[1:] >= 0))[0]
    p1 = None
    if len(j):
        j = j[0]
        d1 = dl[j] - lp[j] / (lp[j + 1] - lp[j]) * (dl[j + 1] - dl[j])
        p1 = (d1, float(np.interp(d1, dl, b.S_num.to_numpy())))
    return mh0, mt_crit, b, p1


def fig_approximation_errors(df, meta, window, out, title, cut):
    """Left: along the M_h line of `cut` nearest the measured M_h (plus the near-critical line at 125.15 GeV), the
    numerical action and the relative errors of the closed forms versus the distance from the stability boundary.
    Right: map of the S_FL error over `window` from `df`."""
    mh0, mt_crit, line, p1 = critical_cut(cut)
    sets = [(line, dict(marker="s", ms=2.5), rf"{cut.attrs.get('label', 'scan')}, $M_h={mh0:g}$ GeV")]
    if os.path.exists(NEAR_CRITICAL):
        nc = pd.read_csv(NEAR_CRITICAL)
        nc = nc[(nc.prescription == "cw1") & np.isclose(nc.Mh, 125.15) & (nc.S_num > 0)].sort_values("delta")
        sets.insert(0, (with_e1(nc), dict(marker="o", ms=3.5), r"near-critical line, $M_h=125.15$ GeV"))
    fig = plt.figure(figsize=(15, 7.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.1], height_ratios=[1, 1.25], hspace=0.07, wspace=0.2)
    a1 = fig.add_subplot(gs[0, 0])
    a2 = fig.add_subplot(gs[1, 0], sharex=a1)
    a3 = fig.add_subplot(gs[:, 1])
    for d, kw, lab in sets:
        a1.loglog(d.delta, d.S_num, color="k", ls="-", lw=0.8, label=lab, **kw)
        for k, (c, _) in ESTIMATES.items():
            e = (d[k] / d.S_num - 1).abs()
            ok = np.isfinite(e) & (e > 0)
            a2.loglog(d.delta[ok], e[ok], color=c, ls="-", lw=0.8, **kw)
    delta_meas = MEASURED["Mt"] - mt_crit
    marks = [(delta_meas, "#1f3cff", f"measured $M_t$ ($S\\approx{np.interp(delta_meas, line.delta, line.S_num):.0f}$)")]
    if p1 is not None:
        marks.append((p1[0], "#a00000", r"$\wp=1$"))
        a1.axhspan(p1[1] * 0.9, p1[1] * 1.1, color=COLORS["unstable"], alpha=0.35, lw=0)
        a1.text(0.06, p1[1] * 1.3, rf"$\wp=1$ needs $S\approx{p1[1]:.0f}$", fontsize=8.5, color="#a00000")
    for x, c, lab in marks:
        for a in (a1, a2):
            a.axvline(x, color=c, lw=1.0, ls="--")
        a1.text(x / 1.1, 1.5e8, lab, color=c, fontsize=9, rotation=90, va="top", ha="right")
    for y, t in ((0.01, "1%"), (0.1, "10%"), (1.0, "100%")):
        a2.axhline(y, color="0.5", lw=0.7, ls=":")
        a2.text(9.5, y * 1.15, t, color="0.35", fontsize=8.5, ha="right")
    for k, (c, lab) in ESTIMATES.items():
        a2.loglog([], [], color=c, lw=1.5, label=lab)
    for a in (a1, a2):
        a.grid(True, which="major", ls="--", alpha=0.4, lw=0.6)
    a1.set_ylim(80, 3e8)
    a1.set_ylabel("numerical bounce action $S$")
    a1.tick_params(labelbottom=False)
    a1.legend(loc="lower left", fontsize=8.5, framealpha=0.95, edgecolor="k")
    a1.set_title(rf"Along $M_h\approx{mh0:g}$ GeV (cw1), versus the distance from the stability boundary", fontsize=11)
    a2.set_xlim(1e-4, 12)
    a2.set_ylim(1e-7, 30)
    a2.set_xlabel(r"$\delta=M_t-M_t^{\rm crit}$ in GeV  ($M_t^{\rm crit}$: absolute-stability boundary)")
    a2.set_ylabel(r"relative error $|S_{\rm est}/S-1|$")
    a2.legend(loc="lower left", fontsize=8.5, framealpha=0.95, edgecolor="k")

    has_s = df.status.isin(BOUNCE) & (df.S_num > 0)
    mh, mt, z = MF.grid(df.assign(_e=(df.S_FL / df.S_num - 1).abs().where(has_s)), "_e")
    a3.set_facecolor("#ececec")
    pm = a3.pcolormesh(mh, mt, np.ma.masked_invalid(z), norm=LogNorm(1e-5, 1), cmap="magma_r", shading="nearest",
                       rasterized=True)
    draw_boundaries(a3, df, lw=0.9)
    draw_measured(a3, star_size=11)
    a3.axvline(mh0, color="#1f3cff", lw=0.8, ls="--")
    cb = fig.colorbar(pm, ax=a3, fraction=0.05, pad=0.02)
    cb.set_label(r"$|S_{\rm FL}/S-1|$ (relative error of the Fubini-Lipatov action)")
    a3.text(0.015, 0.985, "grey: no numerical action (stable, Planck-scale dominated, non-perturbative, action bound,\n"
            r"unresolved); black: $\lambda_{\min}=0$ and $\wp=1$; dashed blue: the cut on the left",
            transform=a3.transAxes, fontsize=8, va="top", bbox=dict(fc="w", ec="0.5", alpha=0.9))
    finish(a3, window[:2], window[2:], title, subtitle(meta, split=True))
    fig.savefig(out)
    plt.close(fig)


def counts(s):
    return {k: int(v) for k, v in s.value_counts().items()}


def crossings(line, col):
    """M_t at which `col` changes sign along a line of constant M_h (linear interpolation between grid points)."""
    l = line[line[col].notna()].sort_values("Mt")
    x, y = l.Mt.values, l[col].values
    i = np.nonzero(np.sign(y[:-1]) != np.sign(y[1:]))[0]
    return [round(float(x[j] - y[j] * (x[j + 1] - x[j]) / (y[j + 1] - y[j])), 4) for j in i]


def phase_summary(full, close, path):
    """The numbers quoted in README section 10.7, written to `path` (JSON)."""
    from scipy.interpolate import griddata
    out = {}
    for name, d in (("full_plane_cw1", full), ("closeup_cw1", close)):
        if d is None:
            continue
        b = d[d.status.isin(BOUNCE)]
        err = (b.S_FL / b.S_num - 1).abs()
        num = b.log10_p_num < 0
        unresolved = d[d.status == "solver-failure"]
        out[name] = {
            "points": len(d), "step_GeV": d.attrs["step"], "status": counts(d.status),
            "solver_failures": int(d.fallback.sum() + d.bound.sum() + len(unresolved)),
            "decided_by_shooting_fallback": int(d.fallback.sum()), "decided_by_action_bound": int(d.bound.sum()),
            "unresolved_Mh_Mt": [[round(float(x), 4), round(float(y), 4)] for x, y in unresolved[["Mh", "Mt"]].values],
            "bounce_points": len(b),
            "class_differs_S_FL_literature_scale": int(((b.log10_p_FL_lit < 0) != num).sum()),
            "class_differs_S_FL_bounce_scale": int(((b.log10_p_FL_bscale < 0) != num).sum()),
            "S_FL_error_above_1pct_only_for_S_from": float(b.S_num[err > 0.01].min()),
            "S_FL_error_above_10pct_only_for_S_from": float(b.S_num[err > 0.1].min()),
            "log10_p_max_where_S_FL_error_above_1pct": float(b.log10_p_num[err > 0.01].max()),
            "S_FL_error_max_for_450_S_650": float(err[(b.S_num > 450) & (b.S_num < 650)].max()),
        }
    d = close if close is not None else full
    if d is not None:
        m, win = MEASURED, 2 * d.attrs["step"]
        near = d[((d.Mh - m["Mh"]).abs() <= win) & ((d.Mt - m["Mt"]).abs() <= win)]
        at = lambda col: float(griddata((near.Mh, near.Mt), near[col], (m["Mh"], m["Mt"]), method="linear"))
        meas = {"scan": d.attrs["label"], "Mh": m["Mh"], "Mt": m["Mt"], "S_num": at("S_num"),
                "log10_p_num": at("log10_p_num"), "log10_Lambda_I_GeV": at("log10_mu_inst")}
        for n in (1, 2):
            a = math.sqrt(SIGMA[n])
            ins = d[((d.Mh - m["Mh"]) / (a * m["dMh"])) ** 2 + ((d.Mt - m["Mt"]) / (a * m["dMt"])) ** 2 <= 1]
            meas[f"{n}sigma_region"] = {"grid_points": len(ins), "status": counts(ins.status),
                                        "S_num_range": [float(ins.S_num.min()), float(ins.S_num.max())],
                                        "log10_p_num_range": [float(ins.log10_p_num.min()),
                                                              float(ins.log10_p_num.max())]}
        out["measured"] = meas
        line = d[(d.Mh - 125.2).abs() < 1e-6]
        out["boundaries_at_Mh_125.2"] = {"scan": d.attrs["label"], "lambda_min_zero_Mt": crossings(line, "lambda_min"),
                                         "p_one_Mt": crossings(line, "log10_p_num")}
    with open(path, "w", newline="\n") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print("summary written to", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", default="results/full_plane_cw1.csv")
    ap.add_argument("--closeup", default="results/closeup_cw1.csv")
    ap.add_argument("--zoom", default="results/controlled_cw1.csv")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--summary", default="results/phase_diagrams_summary.json")
    args = ap.parse_args()
    out = os.path.join(ROOT, args.outdir)
    os.makedirs(out, exist_ok=True)
    full, mfull = load(args.full)
    close, mclose = load(args.closeup)
    zoom, mzoom = load(args.zoom)
    for d, m, lab in ((full, mfull, "full-plane scan"), (close, mclose, "close-up scan"),
                      (zoom, mzoom, "controlled scan")):
        if d is not None:
            d.attrs["step"] = m["step"]
            d.attrs["label"] = lab
    if full is not None or close is not None:
        phase_summary(full, close, os.path.join(ROOT, args.summary))
    fine = close if close is not None else (zoom if zoom is not None else full)
    cut = close if close is not None else zoom      # finest scan of the controlled region, for the error figures
    if cut is not None and full is not None:
        fig_approximation_errors(full, mfull, (0, 250, 0, 250), os.path.join(out, "approximation_errors_full.png"),
                                 "Accuracy of the analytical bounce action: full plane", cut)
    if cut is not None and close is not None:
        fig_approximation_errors(close, mclose, (110, 140, 160, 185),
                                 os.path.join(out, "approximation_errors_closeup.png"),
                                 "Accuracy of the analytical bounce action: close-up", cut)
    if full is not None:
        fig_full(full, mfull, fine, os.path.join(out, "phase_diagram_full.png"))
        fig_full(full, mfull, fine, os.path.join(out, "phase_diagram_full_contours.png"), contours="action")
        fig_full(full, mfull, fine, os.path.join(out, "phase_diagram_full_instability_scale.png"),
                 contours="instability")
        fig_overlay(full, mfull, (0, 250, 0, 250), os.path.join(out, "overlay_full.png"),
                    "Vacuum stability: analytical vs numerical classification", fine=fine)
    if close is not None:
        w = (110, 140, 160, 185)
        fig_window(close, mclose, w, os.path.join(out, "phase_diagram_closeup.png"),
                   "SM vacuum phase diagram: close-up")
        fig_window(close, mclose, w, os.path.join(out, "phase_diagram_closeup_contours.png"),
                   "Close-up with contours of constant action", contours="action",
                   levels=[150, 200, 300, 400, 500, 700, 1000, 2000, 5000, 20000])
        fig_overlay(close, mclose, w, os.path.join(out, "overlay_closeup.png"),
                    "Analytical vs numerical classification: close-up", fine_window=True)
    if close is not None:          # the close-up grid is finer and uses the same matching as the other figures
        zoom, mzoom = close, mclose
    if zoom is not None:
        fig_window(zoom, mzoom, CONTROLLED, os.path.join(out, "phase_diagram_zoom.png"),
                   "SM vacuum phase diagram: controlled region", contours="instability",
                   levels=[9, 10, 11, 12, 13, 14, 16, 18], star_size=10, loc="below")
    print("phase diagrams written to", out)


if __name__ == "__main__":
    main()
