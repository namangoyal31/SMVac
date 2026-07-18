# NEW_OUTLINE.md — Revised Manuscript Structure

## Gu Principle

Every section must answer one question the reader is asking at that point. If a section answers a question no reader would ask, cut it. If it leaves a question unanswered, expand it.

**Reader profile:** A phenomenologist who knows the Buttazzo phase diagram, uses the conformal approximation routinely, and wants to know whether it matters.

---

## Revised Structure

### Title (revised)

*Planck-Scale Rescuing of the Electroweak Vacuum:*  
*Exact Numerical Bounce Solver vs. the Conformal Approximation*

### Author / Affiliation / Date — unchanged

---

### Abstract (12–15 lines)
→ See ABSTRACT_PLAN.md for three candidate versions.

---

### §1 — Introduction (×1 current page)

**Current:** 3.5 pages, covers Higgs discovery, tunneling rate, prior literature, conformal approx, this paper's findings, roadmap.

**Changes:**
- Cut the Higgs discovery paragraph to 2 sentences (it is common knowledge).
- Condense the tunneling rate / threshold explanation (Eqs. 1–3 remain, reduce surrounding prose).
- Keep the conformal approximation equation and its limitation.
- State the two questions: *(i) How accurate is the conformal approximation near the boundary? (ii) What is the quantitative impact of Planck-scale operators?*
- Bullet the four findings as a preview.
- Roadmap sentence stays.

**Target length:** 2 pages (was 3.5).

---

### §2 — Theoretical Framework (×1.5 current pages)

**Current:** §II.A (Effective potential, β-functions) + §II.B (Bounce equation) + §II.C (Validity regime).

**Changes:**
- **§2.1 — RG-Improved Potential with Planck-Scale Operators** (condensed)
  - One equation: V_eff = λ(φ)φ⁴/4 + c₆φ⁶/Λ².
  - One sentence: 3-loop SM β-functions evolved RK4 from M_t to M_Pl (cite Mihaila:2012pz, Luo:2002ti).
  - State NNLO matching thresholds (δλ = -0.00313, δy_t = -0.0587).
  - c₆ scan range: {0, 0.1, 1, 10}.
  - **Remove** the explicit β-function expansion (Eq. 9).
- **§2.2 — Bounce Action and the Conformal Approximation** (condensed)
  - Keep the EOM (Eq. 6) and action integral (Eq. 7).
  - The conformal approximation S = 8π²/3|λ_min|.
  - One paragraph on why it fails when φ_turn → M_Pl.
  - **Remove** the derivation prose for thin-wall limit.
- **§2.3 — Phase Classification** (new, pulled from current §IV.D)
  - Stable: λ > 0 to M_Pl.
  - Metastable: λ < 0, S > S_threshold.
  - Unstable: λ < 0, S < S_threshold.
  - Non-perturbative: couplings blow up.
  - S_threshold = 4 ln(T_U μ_inst) with T_U = 4.3×10¹⁷ s.

**Target length:** 2 pages (was ~3).

---

### §3 — Numerical Methodology (×1 current page)

**Current:** §III.A (Overview) + III.B (RK4) + III.C (Quadrature) + III.D (GSS) + III.E (Parameter scan).

**Changes:**
- **§3.1 — Pipeline Summary** (one paragraph)
  - RGE integration → bounce profile → action evaluation → radius optimization.
  - "Three-loop RK4 at ε = 10⁻¹³, Simpson N = 2048, GSS τ = 10⁻¹³."
  - Total numerical error < 10⁻¹¹ (cite error budget, Table 1).
- **§3.2 — Parameter Space Scan**
  - Grid: 500×500 = 250,000 points, M_t ∈ [155,185], M_h ∈ [115,135].
  - c₆ ∈ {0, 0.1, 1, 10} per point.
  - Dense close-up runs for the SM-relevant region.
- **Move §III.B–III.D text** to Appendix A (pedagogical descriptions of RK4, Simpson, GSS).
- **Move error validation figures** to Appendix B (convergence plots and residual histogram).

**Target length:** 1 page (was 3.5).

---

### §4 — Results (×4 current pages, expanded to 5–6)

**§4.1 — Global Phase Diagram** (was §V.A)
- Two-panel figure side-by-side: analytical (Buttazzo-style) vs. numerical with contours.
- One paragraph: "The global structure is similar, confirming our RGE implementation reproduces Buttazzo:2013uya."
- One paragraph: "The numerical boundary is systematically shifted upward."

**§4.2 — Breakdown of the Conformal Approximation** (was §V.B)
- Scatter plot (S_exact vs S_approx) as the central evidence figure.
- Three claims with supporting evidence:
  1. Global accuracy: 74.9% within 5%.
  2. Systematic failure near boundary: 20–30% mean discrepancy.
  3. Spatial localization: discrepancy concentrated along λ_min = 0 contour (heatmap).
- Sub-figure: discrepancy vs |λ_min| scaling with error bands.

**§4.3 — Phase Boundary Reclassification** (was §V.C + §V.D overlay)
- Primary headline result.
- The strip: 110 points, 0.3–0.5 GeV in M_t.
- Overlay figure (full + close-up).
- "This shift is comparable to the 1σ uncertainty on M_t."
- Explicit statement: the shift has two components — (i) exact solver vs conformal approx at c₆ = 0, (ii) additional shift from c₆ > 0. Provide numbers for each.

**§4.4 — Breakdown of the Semi-Classical Approximation: The Sandwiched Unstable Region** (was §V.E, expanded)
- Three-panel figure (NEW): potential cartoons + virial diagnostic.
- Physical explanation (no longer in Discussion — it is a result).
- Virial ratio S_kin/|S_pot| ≠ 2 as quantitative diagnostic.
- "S_E < 0 signals Fubini ansatz breakdown, not classical roll-down."
- Classification rule: all S ≤ 0 → unstable.

**§4.5 — Universality Across c₆** (was §V.F)
- Two-panel figure (full + zoom).
- Three observations: monotonic shift, ~1.5–2 GeV for c₆ = 10, universal direction.
- "Any UV completion with c₆ > 0 shifts the boundary in the rescuing direction."

**Target length:** 5–6 pages (was ~4).

---

### §5 — Discussion (×1 current page, expanded to 2)

**Current:** §VI.A (Physical implications) + VI.B (Comparison with existing work) + VI.C (Limitations and future directions).

**Changes:**
- **§5.1 — Physical Implications of the Boundary Shift**
  - "A 0.3–0.5 GeV shift is comparable to δM_t ≈ 0.7 GeV → conformal approximation alone cannot determine the SM vacuum's status."
  - "Exact numerical bounce solvers should become standard for precision vacuum stability analyses."

- **§5.2 — Comparison with Previous Work** (NEW table)
  - Explicit comparison table: Buttazzo:2013uya, Branchina:2013jra, Branchina:2014rva, Lalak:2014qua.
  - Columns: RG order, bounce method, Planck ops included, boundary location at M_h = 125 GeV (M_t value), grid resolution, error estimate.

- **§5.3 — Limitations and Future Directions**
  - No one-loop CW corrections (largest missing piece, ~same order as the shift).
  - c₆ as free parameter (no UV completion).
  - Zero-temperature only.
  - Non-perturbative region requires lattice.

**Target length:** 2–2.5 pages (was ~2, but denser).

---

### §6 — Conclusion (×0.5 current page → ×1 page)

**Current:** One paragraph reiterating the abstract.

**Changes:**
- Summary table (4–5 key numbers).
- "The conformal approximation is insufficient within the current 1σ experimental uncertainty on M_t."
- "This work establishes a rigorous numerical baseline for future studies."

---

### Acknowledgments — unchanged.

---

### Appendix A — Numerical Methods Details
- §III.B (RK4) text + adaptive step-size details.
- §III.C (Simpson quadrature) text + compactification mapping.
- §III.D (Golden-Section Search) text.

### Appendix B — Convergence Validation
- Fig. B1: RK4 convergence (was error_rk4_convergence.png).
- Fig. B2: Quadrature convergence (was error_quadrature_convergence.png).
- Fig. B3: GSS plateau (was error_golden_section.png).
- Fig. B4: EOM residual histogram (was error_residual_histogram.png).
- Table B1: Consolidated error budget (warning: not Table 1 in main text — this is a duplicate. Actually: keep Table 1 in main text, move the figures to Appendix B.)

**Correction:** Table 1 stays in main text (§3.1 or §4.1). Only the four convergence *figures* go to Appendix B. The EOM residuals paragraph stays in the main text but is shortened to 2 sentences.

---

## Revised Section Order (Page Budget)

| Section | Current pages | Target pages | Change |
|---------|--------------|--------------|--------|
| Title/Abstract | 1 | 0.5 | Shorter abstract |
| §1 Introduction | 3.5 | 2 | Condense Higgs discovery, roadmap |
| §2 Theory | 3 | 2 | Cut β-function expansion |
| §3 Numerical | 3.5 | 1 | Move details to Appendix |
| §4 Results | 4 | 5.5 | Expand sandwiched region |
| §5 Discussion | 2 | 2.5 | Add comparison table |
| §6 Conclusion | 0.5 | 1 | Add summary table |
| **Total main** | **17.5** | **14.5** | **Tighter narrative** |
| Appendix A | — | 1.5 | Numerical methods |
| Appendix B | — | 1 | Convergence figures |
| **Total** | **17.5** | **17** | ≈ same length, better focus |

## Key Changes from Current Draft

1. **Title shifts from methodology to physics.**
2. **Abstract cuts from 28 lines to 12–15.**
3. **Numerical methodology collapses from 3.5 pages to 1** (details → Appendix A).
4. **Error validation figures move to Appendix B** (only the table stays).
5. **β-function expansion removed** (Eq. 9 cut).
6. **Sandwiched region moves from Results sub-subsection to standalone §4.4** with a new 3-panel figure.
7. **Comparison table added** to Discussion.
8. **Conclusion gets a numbers table.**
9. **Pure SM (c₆ = 0) boundary shift is disentangled** from the c₆ > 0 shift.
10. **Two close-up phase diagrams removed** (redundant with overlay).
