# FIGURE_REVIEW.md — Figure-by-Figure Critique

## Current Inventory (19 figure blocks in the paper)

| # | File | Section | Type | Verdict |
|---|------|---------|------|---------|
| 1 | error_rk4_convergence.png | IV.A | Error validation | **MOVE to Appendix** |
| 2 | error_quadrature_convergence.png | IV.B | Error validation | **MOVE to Appendix** |
| 3 | error_golden_section.png | IV.C | Error validation | **MOVE to Appendix** |
| 4 | error_residual_histogram.png | IV.D | Error validation | **MOVE to Appendix** |
| 5 | analytical_stability_plot.png | V.A | Phase diagram | **KEEP** (Buttazzo-style reference) |
| 6 | numerical_stability_contours.png | V.A | Phase diagram | **KEEP** (headline result) |
| 7 | error_S_scatter.png | V.B | Scatter plot | **KEEP** (key evidence figure) |
| 8 | error_vs_lambda.png | V.B | Scaling plot | **KEEP** (supports scatter) |
| 9 | error_discrepancy_heatmap.png | V.B | Heatmap | **KEEP** (spatial visualization) |
| 10 | error_boundary_reclassification.png | V.C | Phase diagram | **KEEP** (headline result figure) |
| 11 | analytical_stability_closeup.png | V.C | Phase diagram | **REMOVE** (redundant with overlay) |
| 12 | numerical_stability_closeup.png | V.C | Phase diagram | **REMOVE** (redundant with overlay) |
| 13 | overlay_plot.png | V.C | Phase diagram | **KEEP** |
| 14 | overlay_closeup_plot.png | V.C | Phase diagram | **KEEP** |
| 15 | c6_phase_diagram.png | V.E | Phase diagram | **KEEP** |
| 16 | c6_sweep_plot.png | V.E | Phase diagram | **KEEP** |
| — | *(missing)* | V.D | Schematic | **ADD** (potential cartoon for sandwiched region) |

## Recommended Reorganization

### Main Text (6 figures + 1 table)

**Figure 1 — Global Phase Diagrams** (was #5 + #6)
Two-panel side-by-side: (a) Analytical conformal, (b) Exact numerical with constant-action contours. Same color scheme. Same axes. Direct visual comparison.

**Figure 2 — Scatter + Discrepancy** (was #7)
Keep the three-panel S_exact vs S_approx layout. This is the most information-dense figure and should remain prominent.

**Figure 3 — Discrepancy Scaling** (was #8 + #9)
Two-panel: (a) discrepancy vs |λ_min| (error bands), (b) discrepancy heatmap in (M_t, M_h) plane. These tell different but complementary stories.

**Figure 4 — Boundary Reclassification** (was #10 + #14)
Two-panel: (a) full overlay, (b) close-up overlay. The full overlay gives context; the close-up shows the quantitative width.

**Figure 5 — Sandwiched Unstable Region** (NEW)
Three-panel: (a) zoom of phase diagram at M_h ~ 50 GeV showing the sandwiched region, (b) potential V(φ) at three M_t values (one in each region), (c) virial ratio S_kin/|S_pot| as a function of M_t at fixed M_h = 50 GeV. This visually explains the most novel physics.

**Figure 6 — c₆ Universality** (was #15 + #16)
Two-panel: (a) full phase diagram with multiple c₆ boundaries, (b) zoom near SM point. Label explicitly with "c₆ = 0, 0.1, 1, 10".

**Table I — Error Budget** (was Table 1 in §IV.E)
Keep as is. It's clean and informative.

### Appendix Figures (moved from main text)

**Appendix A — Numerical Methods Description**
Move §III.B–III.D text here. Include no figures — the table suffices.

**Appendix B — Convergence Validation** (was #1, #2, #3, #4)
Two-panel page: (top) RK4 + quadrature convergence, (bottom) GSS plateau + residual histogram. These are valuable for specialist readers but break the narrative.

## Specific Plotting Issues

### Color Schemes
- **Stable/Unstable colors:** The current green/yellow/red/black follows Buttazzo convention. This is correct and should not be changed.
- **Discrepancy heatmap:** Switch from the current colormap to a diverging colormap (e.g., RdBu or coolwarm centered at 0) with white at 0% discrepancy. The current scheme implies zero is at one end of the colorbar, which is misleading since the discrepancy is strictly positive.
- **Constant-action contours (#6):** Use a perceptually uniform colormap (e.g., viridis) for the background phase colors, with the contour lines in a contrasting color (white or black with high contrast). Current blue lines on colored background may be hard to read in print.

### Axis and Label Issues
- **All phase diagrams:** Add explicit contour for S = S_threshold(μ_inst) as a dashed line. Currently the boundary is implicit.
- **Scatter plot (#7):** Panel (c) needs the threshold value S_threshold more clearly marked — add a horizontal dashed line.
- **Overlay plots (#13, #14):** The disagreement dots are 8px blue on a 200 dpi figure — at print resolution (600 dpi) they become invisible. Increase marker size or use a small filled circle with a white edge.
- **c₆ sweep (#16):** The legend lists colors but the boundary lines are dashed styles. Use both color + linestyle differentiation for accessibility (print may be grayscale).

### Figures to Replot from Python Scripts

The plotting scripts (`plot_analytical.py`, `plot_numerical.py`, `plot_overlay.py`, `plot_numerical_contours.py`) need the following updates:

1. **Consistent DPI:** Scripts use 200–300 dpi. For PRD/JHEP, export at 600 dpi or as PDF/EPS.
2. **Consistent styling:** `plot_analytical.py` and `plot_numerical.py` use `alpha=0.5` while `plot_numerical_contours.py` uses `alpha=0.3`. Unify.
3. **Serif fonts:** The overlay scripts use `fontsize=12` without serif configuration. Enable `plt.rcParams['font.family'] = 'serif'` with Computer Modern for TeX-matching output.
4. **Tight layout errors:** `plot_overlay.py` calls `fig.tight_layout()` after setting axis limits, which can clip the experimental star marker. Use `constrained_layout` instead or call `tight_layout` before `set_xlim`/`set_ylim`.
5. **Contour bleeding:** `plot_numerical_contours.py` downsamples stable points to 10% to avoid slow triangulation, which causes contours to bleed into the stable region. Use a Delaunay mask or clip the contour domain explicitly.

### Figures to Remove

- **analytical_stability_closeup.png (#11)** and **numerical_stability_closeup.png (#12):** These are superseded by overlay_closeup_plot.png (#14), which shows both boundaries simultaneously. Having three close-ups of the same region is redundant.
- **error_golden_section.png (#3):** The information content (action saturates at τ < 10⁻²) fits in one sentence of text plus the error budget table. The figure shows a plateau — interesting for a methods appendix, extraneous in the main narrative.

### Missing Figure Priority

1. **HIGH: Sandwiched-region schematic** — a 3-panel figure showing potentials at fixed M_h = 50 GeV for M_t in the upper metastable, sandwiched unstable, and lower metastable regions. This figure makes the virial argument visual and is essential for §V.D.
2. **MEDIUM: Pure SM boundary comparison** — a line plot comparing the stability boundary for c₆ = 0 (exact numerical) vs. conformal approximation vs. Buttazzo:2013uya. This disentangles the φ⁶ effect from the exact-solver effect.
3. **LOW: Virial ratio histogram** — distribution of S_kin/|S_pot| across parameter space, showing the breakdown at |λ_min| → 0.
