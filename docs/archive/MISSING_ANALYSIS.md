# Missing Analyses — What the Paper Needs But Does Not Have

**Context:** This document catalogs all analyses that are referenced, implied, or necessary for the paper's claims but are absent from the repository. Each entry is assigned a priority (Critical/High/Medium/Low) for the rewrite.

---

## M1. Convergence Study for the "Exact" Numerical Solver

**Priority: Critical**

**What's missing:** The paper claims "total numerical error < 10⁻¹¹ relative action" but provides:
- No convergence plot (action vs RK4 step size)
- No convergence plot (action vs Simpson N)
- No convergence plot (action vs GSS tolerance)
- No comparison of fixed-step (dt=0.1) vs adaptive RK4
- No comparison of trapezoidal (current) vs Simpson quadrature for V_eff

**What to do:**
1. Run `solver_numerical.cpp` at a single representative point (M_h=125, M_t=173) with dt = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5}
2. Run with Simpson N = {64, 128, 256, 512, 1024, 2048, 4096}
3. Run with GSS tol = {10⁻⁶, 10⁻⁸, 10⁻¹⁰, 10⁻¹², 10⁻¹⁴}
4. Generate a convergence table for the appendix
5. Estimate true error via Richardson extrapolation

---

## M2. Virial Ratio Diagnostic (S_kin / |S_pot|)

**Priority: Critical**

**What's missing:** The paper claims Fubini-Lipatov breakdown is detectable via the virial ratio S_kin/|S_pot| ≠ 2, but:
- S_kin and S_pot are **not exported** in the CSV output (only S_exact = S_kin + S_pot)
- No figure shows this ratio across the parameter space
- No threshold is established for when the ratio deviates from 2

**What to do:**
1. Add `S_kin, S_pot` columns to the solver output
2. Run a targeted scan at fixed M_h = {5, 10, 20, 50, 100} GeV across M_t
3. Plot virial ratio S_kin/|S_pot| vs M_t for each slice
4. Identify where the ratio deviates significantly from 2 (say, >10%)
5. Correlate with negative S_exact and with the stability boundary

---

## M3. c₆ Parameter Sweep

**Priority: Critical**

**What's missing:** The paper's primary physical result (c₆ rescues ~110 points / 0.3–0.5 GeV) has:
- No CSV data with c₆ as a varying parameter
- No run script for the c₆ sweep solver
- No figure showing M_t boundary vs c₆
- No convergence study for the c₆ solver

**What to do:**
1. Wire `archive/solver_c6_sweep.cpp` into the pipeline with a `run_c6_sweep.py` driver
2. Choose a grid: M_h ∈ [115, 135] GeV, c₆ ∈ {0, 0.1, 0.5, 1.0, 2.0, 5.0} (at least 6 values)
3. For each c₆, scan M_t to find the stability boundary
4. Plot ΔM_t_boundary vs c₆ (the headline figure)
5. Count points rescued between c₆=0 and c₆>0
6. Verify "universality" by testing different M_h values

---

## M4. Fractional Error Map (S_exact vs S_approx)

**Priority: High**

**What's missing:** The 20–30% error claim is presented as a generic statement, but:
- No figure shows the fractional error |S_exact - S_approx|/S_exact across the (M_h, M_t) plane
- The claim cannot be verified from any existing figure
- Existing plots show only the phase classification and contour lines

**What to do:**
1. From `plot_overlay.py` or a new script, compute fractional error for each grid point
2. Make a 2D colormap of fractional error (log scale, or linear capped at 100%)
3. Overlay stability boundary contours
4. Quantify: what fraction of the plane has <5% error? 10%? 20%?
5. This is **Fig. 1** candidate — establishes the motivation for the "exact" solver

---

## M5. Bounce Profile Comparison

**Priority: High**

**What's missing:** The paper contrasts "conformal approximation" with "exact numerical" but:
- No figure shows the actual bounce profiles φ(r) from each method
- The Fubini-Lipatov profile is never displayed
- The conformal profile is never displayed
- There is no visual intuition for why they differ

**What to do:**
1. Extract phi(r) from both solvers at a representative point near the boundary
2. Overlay them in a figure (log-log plot)
3. Superimpose the exact potential V_eff(φ) and mark the tunneling endpoints
4. This builds physical intuition and demonstrates the difference visually

---

## M6. Adaptive vs Fixed-Step RK4 Comparison

**Priority: High**

**What's missing:** The code defines `rk4_adaptive_step` but uses fixed dt=0.1:
- No comparison exists between the two
- The paper may not even mention this discrepancy
- A referee will notice the unused adaptive function

**What to do:**
1. Run both methods at 3 representative points
2. Compare S_exact, runtime, number of steps
3. Determine whether adaptive stepping changes results meaningfully
4. Either remove the adaptive function from code, or justify why fixed-step is adequate

---

## M7. Negative Action Diagnostic

**Priority: High**

**What's missing:** 15,879 points have S_exact < 0. The paper attributes this to the "sandwiched region" but:
- Negative S is spread across the full parameter space (M_t ∈ [98, 234], M_h ∈ [0, 244])
- The vast majority of negative S points occur at high M_t (>231.5), not at low M_h
- No diagnostic distinguishes "Fubini breakdown" from "pure numerical inconsistency"
- No figure maps where S_exact < 0 occurs

**What to do:**
1. Create a map of sign(S_exact) color-coded over the (M_h, M_t) plane
2. Check whether |S_exact| is correlated with known issues (λ_R → 0, coupling blow-up)
3. Determine how many negative-S points are at M_h < 50 (the "sandwiched" claim) vs high M_h
4. Decide whether negative S should be treated as Unstable (current) or as a separate category

---

## M8. Analytical Solver Output Check

**Priority: Medium**

**What's missing:** `results/analytical_data.csv` has only 4 columns (Mt, Mh_calc, Stability, S_approx), no S_exact:
- The overlay script compares M_h between two CSVs with different schemas
- Point counts differ: numerical has 1,737,982 pts, analytical has 941,176 pts
- Matching requires rounding and cKDTree — this is fragile
- No validation that the two solvers sample the same (M_h, M_t) points

**What to do:**
1. Standardize CSV output format between solvers
2. Validate that both solvers run on the same grid
3. Check for point count mismatch source (941k vs 1.7M)

---

## M9. Contour Plot Accuracy

**Priority: Medium**

**What's missing:** `plot_numerical_contours.py` downsamples stable points to 10%:
- This was apparently needed for render time
- But it may cause contour bleeding or misplaced boundaries
- The 10% downsampling factor is not documented in the paper

**What to do:**
1. Check if removing the downsampling changes contour locations
2. Document the downsampling in the figure caption or methods section
3. Consider using a regular grid (imshow/pcolormesh) instead of scatter+contour

---

## M10. NNLO Matching Validation

**Priority: Medium**

**What's missing:** `get_nnlo_matching` uses hardcoded threshold corrections:
- δλ = -0.00313, δyt = -0.0587
- No comparison with `mr` (RGE evolution tool) or with Buttazzo:2013uya
- No validation that these threshold corrections are correct for the solver's RG scale choices

**What to do:**
1. Compare output at M_h=125.1, M_t=173.3 with Buttazzo:2013uya Table I or Fig 1
2. Verify that the stability boundary matches the literature (~M_t = 171–173 GeV for M_h=125)
3. Document the matching scheme and its uncertainties

---

## M11. Experimental Point Overlay with Uncertainty

**Priority: Low**

**What's missing:**
- No error ellipses (1σ, 2σ) for (M_h, M_t) around (125.1, 173.3)
- No discussion of how δM_t ≈ 0.7 GeV affects the stability classification
- The figure shows only a star marker at the central value

**What to do:**
1. Add 1σ and 2σ contours from PDG values
2. Add a short discussion: "At 1σ, the SM vacuum is metastable; at 2σ, stability/instability cannot be excluded"
3. Reference the PDG 2024 values

---

## Priority Summary

| ID | Analysis | Priority | Effort Estimate |
|----|----------|----------|-----------------|
| M1 | Convergence study | **Critical** | 2-3 days |
| M2 | Virial ratio diagnostic | **Critical** | 3-5 days |
| M3 | c₆ parameter sweep | **Critical** | 5-10 days |
| M4 | Fractional error map | **High** | 1 day |
| M5 | Bounce profile comparison | **High** | 1-2 days |
| M6 | Adaptive vs fixed-step RK4 | **High** | 1 day |
| M7 | Negative action diagnostic | **High** | 2-3 days |
| M8 | Analytical solver output check | **Medium** | 0.5 day |
| M9 | Contour plot accuracy | **Medium** | 0.5 day |
| M10 | NNLO matching validation | **Medium** | 1-2 days |
| M11 | Experimental uncertainty overlay | **Low** | 0.5 day |

**Total estimated effort: ~17-29 days** of continuous computation and analysis. This is realistic for a major revision cycle (6-8 weeks).
