# Claim Audit (Skeptical PRD Referee)

**Auditor:** LLM-assisted systematic review  
**Data sources:** `results/numerical_data.csv` (1,737,982 pts), `results/analytical_data.csv` (941,176 pts), `solver_numerical.cpp`, `solver_analytical.cpp`, `paper/main.tex`, `prompt_for_claude.md`  
**Date:** 2026-06-16

---

## Methodology

Each claim is rated:

- ✅ **Supported** — evidence in data/code is consistent with claim
- ⚠️ **Partially supported** — claim is technically true but overstated or incomplete
- ❌ **Not supported** — data contradicts claim or claim cannot be verified from available materials

---

## Claim 1 — "The conformal approximation fails near the metastability boundary, with 20–30% error in the bounce action"

**Rating: ⚠️ PARTIALLY SUPPORTED**

**Evidence:**
- At M_h~125, M_t~173 (experimental point): S_exact ≈ 2433 vs S_approx ≈ 2374 — **error is 2.4%**, far below 20–30%
- Across all 410,625 computed points: only **74.6% match within 1%** — so 25.4% disagree by >1%
- The disagreements concentrate in two regimes: (1) near the 2↔3 boundary where S is close to the threshold, and (2) at very high M_t where S_exact becomes negative

**Referee concern:** The paper states "20–30%" as a general statement, but most of the (M_h, M_t) plane shows much smaller discrepancies. The 20–30% figure applies only to a narrow band around the metastability boundary. This needs to be quantified with a fractional-error map, not described as a generic property. The claim as written in the abstract is misleading.

---

## Claim 2 — "The c₆ φ⁶/M_Pl² operator rescues ~110 points, corresponding to 0.3–0.5 GeV in M_t"

**Rating: ❌ CANNOT BE VERIFIED**

**Evidence:**
- The c₆ sweep solver (`archive/solver_c6_sweep.cpp`) is a **separate executable** not wired into the main pipeline
- No c₆-column CSV data exists in `results/`
- The c₆ solver output format is 4-column (M_t, M_h_calc, c_6, S_exact) — does not match production CSV format
- No `run_c6_sweep.py` script exists
- The reported numbers (110 points rescued, 0.3–0.5 GeV shift) appear in the LaTeX source as hardcoded values with no reference to a supporting analysis script

**Referee concern:** The c₆ claim is one of the paper's two primary physical results, but it is entirely unverifiable from the available data and code pipeline. Even the supporting solver compiles separately and has no run script. This would require significant additional work during review — either a dedicated reproducibility run or a complete re-implementation.

---

## Claim 3 — "The Fubini-Lipatov ansatz breaks down at low Higgs masses, producing a 'sandwiched' unstable region where S_kin/|S_pot| ≠ 2"

**Rating: ❌ NOT SUPPORTED BY DATA**

**Evidence:**
- **No sandwiched (2→3→2) pattern exists** in the CSV data at any M_h slice:
  - M_h ~ 50 GeV: pattern is Metastable (2) → Unstable (3) with rapid 2↔3 oscillations at the boundary, but **no return from 3→2**
  - M_h ~ 70 GeV: similar 2→3 with isolated negative-S outliers at low M_t
  - M_h ~ 125 GeV: entirely Unstable (3) for M_t > 171 GeV
- **Negative S_exact** (15,879 points) is found **across the entire parameter space** (M_t ∈ [97.8, 233.5], M_h ∈ [0, 244]), not localized to a sandwiched region at low M_h
- The most negative S_exact occurs at **high** M_t and M_h (e.g., M_t=229, M_h=236.5: S = -8.13×10⁸), not at low M_h
- **S_kin and S_pot are not exported** in CSV output. The virial ratio S_kin/|S_pot| ≠ 2 cannot be verified from any available data
- Negative action arises from the **inconsistency** of using a constant-λ kinetic term formula (16π²/3|λ_R|) with a running-λ potential integral — not from actual barrier disappearance

**Referee concern:** This claim appears to be the paper's secondary physical discovery, but it is entirely unreproducible. The virial-ratio diagnostic is not in the data. The claimed "sandwiched" region is not visible in the actual stability classifications. The negative action is a methodological artifact, not a physical discovery. This section would need to be substantially revised or removed.

---

## Claim 4 — "We use an exact numerical bounce action solver with RK4 tolerance 10⁻¹³ and total numerical error < 10⁻¹¹ relative action"

**Rating: ❌ MISLEADING/MISCHARACTERIZED**

**Evidence:**
- **The solver is NOT an exact numerical bounce solver** — it uses the Fubini-Lipatov profile ansatz:
  - `phi(r) = sqrt(2/|lambda_R|) * 2R / (r^2 + R^2)` (line ~207)
  - `S_kin = 16*pi^2 / (3*|lambda_R|)` (line ~230, closed form, assumes constant λ)
  - `S_pot` is the only numerically integrated term, using trapezoidal quadrature (line ~245)
- The solver is a **semi-analytical hybrid**, not a full ODE bounce solution
- **rk4_adaptive_step is defined but NEVER CALLED** (lines ~162-196)
- The main RK4 integration uses **fixed step size dt=0.1** (line ~283)
- Convergence tolerance 10⁻¹³ applies only to the **GSS radius optimization**, not to the RK4
- Total numerical error < 10⁻¹¹ is asserted without any convergence study shown in the code

**Referee concern:** The paper frames this as a full numerical calculation throughout the abstract and introduction. The word "exact" appears multiple times. This is a mischaracterization — the method is an ansatz-based semi-analytical evaluation. A referee familiar with bounce calculations will immediately identify the Fubini-Lipatov profile and the closed-form kinetic term and realize this is NOT the advertised method. This is the single most damaging issue for the paper.

---

## Claim 5 — "The effect of the c₆ operator is universal — any positive c₆ shifts the metastability boundary upward"

**Rating: ⚠️ LOGICALLY SOUND BUT UNVERIFIED**

**Evidence:**
- This is a statement about the RG running of the Higgs quartic: a positive c₆ term adds a positive contribution to the β-function of λ at high energies, which raises λ(M_Pl)
- This is theoretically well-motivated and consistent with known physics
- **No c₆ scan data exists in the repository** to demonstrate the claimed universality
- The c₆ sweep solver (`archive/solver_c6_sweep.cpp`) exists but was never run or integrated into the pipeline

**Referee concern:** The statement is physically plausible but entirely unsubstantiated in the manuscript. No figure shows the claimed universal effect. No grid scan with varying c₆ is presented. The paper would benefit from removing or substantially reducing this claim, or actually performing the analysis.

---

## Claim 6 — "500×500 grid (250,000 points) covering M_t ∈ [155, 185], M_h ∈ [115, 135]"

**Rating: ⚠️ INACCURATE BUT CONSERVATIVE**

**Evidence:**
- The CSV contains **1,737,982 total points**, encompassing a much wider range (M_t ∈ [~94, 236], M_h ∈ [~0, 246])
- Of these, 410,625 have computed S_exact (non-sentinel)
- The 250k figure in the paper is a subset — the actual computation is **7× larger**

**Referee concern:** This is a minor inaccuracy that actually undersells the work. However, the CSV was apparently accumulated from multiple runs, so the effective grid for any single analysis may indeed be ~250k. This should be clarified.

---

## Claim 7 — "The Standard Model vacuum is metastable: life-time exceeds the age of the universe"

**Rating: ✅ SUPPORTED**

**Evidence:**
- At M_h≈125.1, M_t≈173.3: S_exact ≈ 1860–2433 (Stab=2, Metastable)
- The decay threshold S_threshold ≈ 4 ln(T_U μ_inst) ≈ 250–340 for typical μ_inst
- S_exact ≫ S_threshold with comfortable margin
- At higher M_t (~171): S drops below threshold, becoming Unstable (Stab=3)
- The stability boundary lies at M_t ≈ 171–173 for M_h≈125, consistent with Buttazzo:2013uya

**Referee assessment:** This is the most robust result in the paper. It confirms the well-known vacuum metastability result and adds a more precise boundary determination. This claim is fully supported.

---

## Summary Table

| # | Claim | Rating | Severity |
|---|-------|--------|----------|
| 1 | Conformal approx 20–30% error near boundary | ⚠️ Partial | Medium |
| 2 | c₆ rescues ~110 points / 0.3–0.5 GeV | ❌ Cannot verify | **High** |
| 3 | Fubini breakdown → sandwiched unstable region | ❌ Contradicted | **High** |
| 4 | Exact numerical solver, tol 10⁻¹³ | ❌ Mischaracterized | **Critical** |
| 5 | c₆ effect universal | ⚠️ Unverified | Medium |
| 6 | 250k grid points | ⚠️ Conservative | Low |
| 7 | SM vacuum metastable | ✅ Supported | — |

---

## Verdict

**Two claims are supported by data:** the metastability of the SM vacuum (claim 7) and the superiority of the numerical approach over the conformal approximation (claim 1, partially).

**Two claims are fatal as written:** the solver methodology is mischaracterized (claim 4) and the sandwiched unstable region is not visible in the data (claim 3).

**One claim is unverifiable from available materials:** the c₆ rescaling (claim 2).

The paper needs major revisions before it can pass peer review:
1. The introduction must honestly describe the Fubini-Lipatov ansatz method
2. The sandwiched region claims must be either removed or substantiated with the virial ratio
3. The c₆ results must be either removed or actual sweep data must be produced
4. The 20–30% error claim must be qualified with a fractional-error map
