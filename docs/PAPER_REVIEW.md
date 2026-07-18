# PAPER_REVIEW.md — Editorial & Scientific Assessment

## Classification

This is a **combination paper** with four interleaved contributions:

| Component | Weight | Novelty |
|-----------|--------|---------|
| **Phase-boundary correction** (0.3–0.5 GeV shift in M_t) | Primary | High — first exact quantification |
| **Fubini/virial-breakdown discovery** (sandwiched unstable region) | Secondary | High — new physical mechanism |
| **Precision numerical methods** (sub-10⁻¹¹ error pipeline) | Enabling | Moderate — known methods, first application at this scale |
| **Vacuum topology** (constant-action contours) | Supporting | Moderate — novel visualization |

## Headline Result

**"The conformal approximation cannot resolve whether the SM vacuum is metastable or unstable within 1σ of M_t — exact numerical solvers are strictly necessary."**

This is the strongest claim because it speaks directly to the experimental community (ATLAS/CMS top mass program) and the theoretical community (vacuum stability literature) simultaneously.

## What Needs Attention

### Sections That Belong in Appendices

- **§III.B–III.D (RK4, Quadrature, GSS descriptions):** Move to Appendix A. The error budget table (Table I) belongs in the main text; the pedagogical descriptions of Simpson's rule and Golden-Section search do not.
- **§IV.A–IV.C (convergence plots):** Move to Appendix B. Keep only the error budget table in the main text. RK4/quadrature/GSS convergence figures are specialist validation that interrupts the narrative flow.

### Sections That Should Be Expanded

- **§IV.D (EOM Residuals):** Currently one paragraph + histogram. This is the strongest direct validation of the solver. Add a brief discussion of where residual outliers occur (near stability boundary) and why they are harmless.
- **§V.B.3 (Sandwiched unstable region):** This is the most novel physical result in the paper but gets only 1.5 pages. Needs a schematic figure (cartoon of the potential vs. M_t) and a clearer statement distinguishing Fubini-ansatz breakdown from true instability.
- **§VI (Discussion):** The comparison with existing literature is thin. Add a dedicated table comparing this work's numbers (boundary location, shift width, S_E values) with Buttazzo:2013uya, Branchina:2013jra, etc.

### Sections That Should Be Cut or Condensed

- **§II.A (RG-improved effective potential):** The three-loop β-function formula (Eq. 9) takes a half page but is never referenced again. State that 3-loop SM β-functions are used and cite Mihaila:2012pz. Move the explicit expansion to an appendix if needed.
- **§II.B (Bounce equation):** Standard textbook material. Keep the EOM and action integral, cut the derivation prose.

### Missing Content

- **No explicit comparison table** of boundary positions (M_t at fixed M_h) between this work and Buttazzo:2013uya.
- **The c₆ = 0 (pure SM) exact numerical result** is never shown independently of the analytical comparison. Is the boundary shift purely from the φ⁶ operator, or does the exact solver also shift the pure-SM boundary relative to the conformal approximation? This must be disentangled.
- **No error bars on the boundary location.** The 110-point / 0.3–0.5 GeV strip width is stated without uncertainty. A bootstrap or jackknife estimate is needed.

### Recommendations

1. **Reframe the title** to emphasize the boundary shift over the methodology: *"Exact Numerical Bounce Solver Reveals Planck-Scale Rescuing of the Electroweak Vacuum"* or similar.
2. **Shorten the abstract.** The current abstract (28 lines) reads like a summary. Cut to 12 lines: state the problem, the method (one sentence), the headline number (0.3–0.5 GeV shift, 110 points), and the implication.
3. **Merge §V (Results) subsections.** The current 5 subsections + 2 sub-subsections creates fragmentation. Restructure as:
   - V.A Global Phase Diagram (analytical + numerical side-by-side)
   - V.B Breakdown of the Conformal Approximation (scatter + λ_min correlation + heatmap)
   - V.C Boundary Reclassification (overlay + closeup + strip characterization)
   - V.D Sandwiched Unstable Region (virial diagnostic + potential cartoon)
   - V.E c₆ Universality
4. **Add a Conclusions table** summarizing the 4-5 key numbers (boundary shift in GeV, number of reclassified points, relative discrepancy range, c6 sensitivity).
