# FIGURE_PLAN.md — Final Figure Order and Specifications

## Master Table

| Order | Figure ID | Content | Panel layout | Status |
|-------|-----------|---------|--------------|--------|
| 1 | Global phase diagrams | Analytical + numerical with contours | 2-panel (a)(b) side-by-side | **Replot** |
| 2 | S_exact vs S_approx scatter | Three-panel comparison | 3-panel (a)(b)(c) | Keep, **reformat** |
| 3 | Discrepancy scaling | (a) Discrepancy vs \|λ_min\| (b) Heatmap | 2-panel (a)(b) | **Replot heatmap** |
| 4 | Boundary reclassification | (a) Full overlay (b) Close-up overlay | 2-panel (a)(b) | **Replot markers** |
| 5 | Sandwiched unstable region | (a) Phase zoom (b) Potential cartoon (c) Virial ratio | 3-panel (a)(b)(c) | **NEW (must create)** |
| 6 | c₆ universality | (a) Full boundaries (b) SM-region zoom | 2-panel (a)(b) | **Replot** |
| — | Table I | Error budget | — | Keep as-is |
| B1 | RK4 convergence | log-log error scaling | single panel | Move to Appendix |
| B2 | Quadrature convergence | log-log error scaling | single panel | Move to Appendix |
| B3 | GSS plateau | S_min vs tolerance | single panel | Move to Appendix |
| B4 | EOM residual histogram | Residual distribution | single panel | Move to Appendix |

---

## Figure 1 — Global Phase Diagrams

**File:** `figures/global_phase_diagrams.png` (new combined plot)

**Layout:** 2 panels, side-by-side, same aspect ratio (~10×4 total).

**Panel (a):** `analytical_stability_plot.png` — as currently generated.
**Panel (b):** `numerical_stability_contours.png` — as currently generated.

**Changes needed in plotting scripts:**
- Both panels must share identical axis ranges: x = [0, 250], y = [0, 250].
- Both panels must use the same point size and alpha for scatter points.
- Add a horizontal dashed line in panel (b) at M_t = 173.3 GeV and a vertical dashed line at M_h = 125.1 GeV, labeling them "M_t^{exp}" and "M_h^{exp}".
- The experimental star must appear in both panels.
- Panel (b) contour labels: use black text with white halo for readability.
- Supertitle: "Figure 1. Global vacuum stability phase diagrams."
- Panel labels: "(a) Conformal approximation" and "(b) Exact numerical (c₆ = 1)".

**Remove** the existing standalone `analytical_stability_plot.png` and `numerical_stability_contours.png` from the paper (they are now sub-panels).

---

## Figure 2 — Bounce Action Comparison

**File:** `figures/S_scatter.png` (exact copy of `error_S_scatter.png`)

**Layout:** 3 panels (a)(b)(c) as currently designed — this figure is excellent and should be preserved nearly verbatim.

**Changes needed:**
- Panel (c): Add a horizontal dashed line at S = S_threshold ≈ 400 with label "S_threshold ≈ 400".
- Panel (c): Increase the marker size of the rescued points (green box) so they are visible at print resolution.
- Add the 74.9% and 86.8% statistics as text annotations in panel (a), not just in the caption.

---

## Figure 3 — Discrepancy Scaling

**File:** `figures/discrepancy_scaling.png` (new combined plot)

**Layout:** 2 panels (a)(b), vertically stacked or side-by-side.

**Panel (a):** `error_vs_lambda.png` — as currently generated.
- Add a horizontal band at ±5% to visually separate "acceptable" from "failed" regimes.
- Annotate the "Boundary region" where |λ_min| < 10⁻³ with a vertical shaded band.

**Panel (b):** `error_discrepancy_heatmap.png` — **replot with diverging colormap**.
- Current issue: the colormap has zero-discrepancy at one end, which is misleading (the discrepancy is strictly positive in this problem, but the colorbar suggests zero is a boundary).
- Fix: switch to `coolwarm` (or `RdBu`) centered at 0, with 0 mapped to white. Since all values are positive, the entire heatmap will appear in the warm half — but crucially the colorbar will correctly show 0 at the midpoint, not at the minimum.
- Add M_t^{exp} and M_h^{exp} dashed guidelines.

---

## Figure 4 — Boundary Reclassification

**File:** `figures/boundary_reclassification.png` (new combined plot)

**Layout:** 2 panels (a)(b), vertically stacked.

**Panel (a):** Full overlay. `overlay_plot.png` as currently generated.
**Panel (b):** Close-up overlay. `overlay_closeup_plot.png` as currently generated.

**Changes needed:**
- **CRITICAL:** The blue disagreement dots in `overlay_plot.py` are 8px at 200 dpi → at 600 dpi they will be ~2–3 px, essentially invisible. Increase to `s_disagree=12` at minimum, or use `s_disagree=20` with `edgecolor='white'` and `linewidth=0.3` for print visibility.
- In panel (b), add a text annotation "∼0.3–0.5 GeV" with a horizontal brace or arrow spanning the width of the blue strip.
- Convert both panels to 600 dpi output.
- Disable `rasterized=True` for the disagreement dots (rasterization causes artifacts at print resolution). Use vector output for these.
- Add the experimental 1σ ellipse for M_t (vertical band at M_t = 173.3 ± 0.7) in panel (b) as a shaded region, to visually show the overlap between the boundary shift and the experimental uncertainty.

**See also:** Remove the standalone close-up figures `analytical_stability_closeup.png` and `numerical_stability_closeup.png` — they are fully superseded by this figure.

---

## Figure 5 — Sandwiched Unstable Region (NEW)

**File:** `figures/sandwiched_region.png`

**Layout:** 3 panels (a)(b)(c), 1 row × 3 columns or 2+1 layout.

This is the only figure that does not currently exist. It requires generating new data/plots.

**Panel (a) — Phase diagram zoom at low M_h:**
- Zoom of the numerical phase diagram: M_h ∈ [0, 120] GeV, M_t ∈ [80, 200] GeV.
- Highlight the "sandwiched" region where the classification is S = 3 (unstable) sandwiched between S = 2 (metastable) bands.
- Overlay the virial ratio as a color dimension on the scatter points (or as a contour).

**Panel (b) — Potential V(φ) at three M_t values (fixed M_h = 50 GeV):**
- Three curves: M_t in the upper metastable band, M_t in the sandwiched unstable band, M_t in the lower metastable band.
- x-axis: φ (linear, 0–10⁴ GeV or log scale as appropriate).
- y-axis: V(φ) (showing barrier, well, and false vacuum).
- Annotate barrier height and width for each curve.
- Use a linestyle+color scheme to distinguish the three curves.
- Add an inset zoom of the false-vacuum region around φ ∼ v.

**Panel (c) — Virial ratio S_kin/|S_pot| vs M_t at fixed M_h = 50 GeV:**
- x-axis: M_t, y-axis: S_kin/|S_pot|.
- Horizontal dashed line at S_kin/|S_pot| = 2 (virial constraint for pure quartic).
- Shaded bands for the three stability regions.
- Annotate the "Fubini breakdown" regime where the ratio deviates from 2 significantly.
- On the right y-axis, optionally overlay the computed S_exact value (log scale) to show the correlation.

**Data needed:**
- The virial ratio S_kin/|S_pot| is already computed in the solver output (internal variables used in action evaluation). The numerical CSV files contain S_exact, S_kin, S_pot — check if S_kin and S_pot are exported. If not, a small modification to the C++ solver (adding S_kin and S_pot to the CSV output) is needed. This is a high-priority task.

---

## Figure 6 — c₆ Universality

**File:** `figures/c6_universality.png` (new combined plot)

**Layout:** 2 panels (a)(b).

**Panel (a):** `c6_phase_diagram.png` — as currently generated.
**Panel (b):** `c6_sweep_plot.png` — as currently generated.

**Changes needed:**
- In both panels, the legend must differentiate by **both color and linestyle** for grayscale accessibility. Use solid/dashed/dotted/dash-dot for c₆ = 0/0.1/1/10 respectively.
- In panel (a), the c₆ = 0 (pure SM) boundary is shown in blue — clearly label it "c₆ = 0 (pure SM)".
- Add a reference: "Buttazzo:2013uya bound" as a dotted gray line on both panels if the c₆ = 0 boundary differs from the literature.
- Panel (b): The experimental star should be the most prominent element (increase zorder, use white fill with thick black edge).
- Add text annotation showing the arrow of increasing c₆ with label "Increasing c₆" in panel (a).

---

## Table I — Error Budget

Same as current Table 1 in §IV.E. Keep format, wording, and placement. It is the only table in the main text and serves as the compact replacement for the four convergence figures.

**Requested addition:** Add a row showing NNLO matching uncertainties (estimated from varying δλ and δy_t within their ranges) to show that parametric uncertainty dominates numerical uncertainty.

---

## Appendix B — Convergence Validation (4 panels)

**File:** `figures/appendix_convergence.png` (new combined plot)

**Layout:** 2×2 grid.

| Top-left | Top-right |
|----------|-----------|
| `error_rk4_convergence.png` | `error_quadrature_convergence.png` |

| Bottom-left | Bottom-right |
|-------------|--------------|
| `error_golden_section.png` | `error_residual_histogram.png` |

All four panels should share the same styling (font, axis label size, grid style).

**Changes needed:**
- Standardize all four plots to use `plt.rcParams['font.family'] = 'serif'` with Computer Modern.
- Use consistent `figsize` per panel.
- The GSS plateau figure (`error_golden_section.png`) currently has no y-axis label. Add: "S_min (relative change)".
- The residual histogram (`error_residual_histogram.png`) x-axis label should read "EOM residual (integrated)".

---

## Color and Style Specification

| Element | Specification |
|---------|--------------|
| Stable | `#78f078` (green) — same as Buttazzo |
| Metastable | `#ffff66` (yellow) — same as Buttazzo |
| Unstable | `#ff7878` (red) — same as Buttazzo |
| Non-perturbative | `#505050` (gray) — same as Buttazzo |
| Experimental star | Blue (`#1f77b4`) with black edge, marker `*`, size 200–300 |
| Disagreement dots | Blue (`#1f77b4`) with white edge, s=15–20, edgewidth=0.5 |
| Contour lines | Black or white (whichever contrasts with background), linewidth=1.5 |
| Contour labels | White text with black halo (`path_effects=[Stroke(fill_color='black')]`) |
| Font | Computer Modern serif (via `rcParams['font.family'] = 'serif'`) |
| DPI | 600 for main figures, 300 for appendix |
| Format | PDF for submission, PNG for preview |

## Figure Order in the Manuscript

1. **Figure 1** — Global phase diagrams — §4.1
2. **Figure 2** — S_exact vs S_approx scatter — §4.2
3. **Figure 3** — Discrepancy scaling and heatmap — §4.2
4. **Figure 4** — Boundary reclassification — §4.3
5. **Figure 5** — Sandwiched unstable region — §4.4 (NEW)
6. **Figure 6** — c₆ universality — §4.5
7. **Table I** — Error budget — §3.1 (or Appendix B reference)

**Appendix B Figures:** B1–B4, referenced as "Fig. B1" through "Fig. B4" in the Appendix.

## Script Modification Requirements

| Script | Changes needed |
|--------|---------------|
| `plot_analytical.py` | Add serif font config; add M_t^{exp}, M_h^{exp} guidelines; 600 dpi; shared axis ranges |
| `plot_numerical.py` | Same changes as analytical |
| `plot_numerical_contours.py` | Fix contour bleeding (Delaunay mask); unify alpha=0.5; serif fonts; 600 dpi |
| `plot_overlay.py` | Increase markersize to 15–20; remove rasterization; add 1σ M_t ellipse; add width annotation; 600 dpi |
| — (new script) | `plot_sandwiched_region.py` — generate Figure 5 from solver virial data |
