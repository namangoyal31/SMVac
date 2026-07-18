# Claim Impact Assessment

**Purpose:** Map each paper claim to its required revision, severity, and concrete
rewrite recommendation, based on the code-level verification in SOLVER_VERIFICATION.md.

**Role:** Hostile PRD referee issuing a mandatory revision report.

---

## Impact Classification

- **FATAL** — Claim is false and central to the paper's contribution. Paper cannot be published without removing or fundamentally rewriting.
- **MAJOR** — Claim is misleading or overstated. Paper can be published after correction and qualification.
- **MINOR** — Claim is technically inaccurate but does not affect the scientific conclusions.

---

## Impact 1: "Exact numerical bounce solver" (Abstract, Intro, §III)

**Severity: FATAL**

**Where it appears:**
- Abstract line 31: "exact numerical framework"
- Abstract line 38: "exact bounce solver"
- Intro line 11: "fully exact numerical pipeline"
- §III header: "Numerical Methodology" (unqualified)
- §III.A: "exact bounce action computation"
- §III.C caption: "exact numerical bounce solver"

**What the code shows:**
The solver uses a Fubini-Lipatov ansatz with closed-form kinetic term.
There is no ODE solve. The method does **not** converge to the true bounce
solution in any limit (see SOLVER_VERIFICATION.md §5).

**Required revision:**
- Remove every instance of "exact" applied to the bounce solver.
- Replace with: *"semi-analytical bounce action computation based on the
  Fubini-Lipatov ansatz with numerically integrated running-coupling corrections"*
- Add a methodology subsection ("Bounce Profile Ansatz") that explicitly states
  the Fubini-Lipatov profile, its closed-form kinetic term, and that only the
  radius and potential integral are evaluated numerically.
- All subsequent references ("exact result", "exact solver", "exact diagram")
  must be replaced throughout.

---

## Impact 2: "Numerical errors below $10^{-11}$ in relative action" (Abstract, §IV)

**Severity: FATAL**

**Where it appears:**
- Abstract line 37: "we achieve numerical errors below $10^{-11}$ in relative action"
- §IV.E, Table I: "Total numerical error $\lesssim 10^{-11}$"
- Captions of Figs. 3–6

**What the code shows:**
The error budget (Table I) includes RK4 truncation, Simpson quadrature, and
GSS optimization errors, but **omits the dominant systematic error**: the
Fubini-Lipatov ansatz approximation itself. Appendix D admits this error can be
$\mathcal{O}(1)$ or divergent (negative action). The $10^{-11}$ figure is the
error of evaluating a **wrong profile accurately**.

**Required revision:**
- Table I must include a row: *"Profile ansatz systematic error (Fubini-Lipatov)"*
  with an estimated magnitude. The current data show this error reaches 25%+ near
  the boundary and produces completely wrong signs (negative action) in $\sim 4\%$
  of computed points.
- The abstract's "$10^{-11}$" claim must be qualified: *"...numerical errors in
  the RGE integration, quadrature, and radius optimization below $10^{-11}$ in
  relative action, with additional systematic uncertainty from the Fubini-Lipatov
  profile ansatz that can reach $\mathcal{O}(10\%)$ near the boundary"*
- Remove or relocate the "$10^{-11}$" claim from the abstract; it is misleading
  as a headline.

---

## Impact 3: "Adaptive step sizing" for the RK4 integrator (§III.A)

**Severity: MAJOR**

**Where it appears:**
- §III.A lines 284–294: "adaptive step-size controller ... automatically refining
  the grid near regions of rapid coupling evolution"

**What the code shows:**
`solver_numerical.cpp:282–304` defines `rk4_adaptive_step` with tolerance $10^{-13}$,
but the function is **never called**. The main loop (`solver_numerical.cpp:449`) calls
`rk4_single_step` with **fixed dt = 0.1**. The step-size never changes.

**Required revision:**
- Replace "adaptive step sizing" with "fixed-step RK4 with $\Delta t = 0.1$"
- Provide a justification for why fixed-step is adequate (e.g., the RGE functions
  are smooth enough that adaptive stepping is unnecessary)
- Either remove the `rk4_adaptive_step` function from the code, or add a note
  explaining why it was not used

---

## Impact 4: "Equation-of-motion residuals confirm bounce profiles satisfy the ODE" (§IV.D)

**Severity: MAJOR**

**Where it appears:**
- §IV.D lines 452–454: "the numerical bounce profiles tightly satisfy their
  defining ODE throughout the parameter space"
- Figure 6 caption and histogram description

**What the code shows:**
The "EOM residual" is computed by plugging the Fubini-Lipatov profile into the
bounce ODE and measuring the mismatch. A non-zero residual does **not** indicate
numerical error — it indicates the **profile is not a solution**. The histogram
peak near zero means the Fubini-Lipatov profile approximately satisfies the ODE
**when $\lambda$ runs slowly**, which is expected from its derivation. The large
tail corresponds to the boundary region where the ansatz fails — which is the
paper's own claimed discovery region.

**Required revision:**
- Reframe: *"We quantify the accuracy of the Fubini-Lipatov ansatz by evaluating
  its residual in the bounce equation of motion"*
- Remove the language implying this validates the numerics
- Add a caveat: large residuals occur precisely where the paper claims new physics
- Consider moving this figure to an appendix

---

## Impact 5: "The conformal approximation has 20–30% error near the boundary" (§II.B, §V)

**Severity: MAJOR**

**Where it appears:**
- Abstract lines 32–34: "this approximation systematically fails near the
  metastability boundary"
- §II.B lines 257–258: "catastrophically underestimates $S_E$"
- §V text and Figure 3 caption: 20–30% error claim
- Line 142: "74.9% of points agree within 1%"

**What the code/data show:**
- At the experimental point (M_h=125, M_t=173.3): S_approx=2374, S_exact=2433,
  **error = 2.4%**, not 20–30%
- 74.6% of computed points match within 1% — consistent with the paper's
  74.9% claim
- The 20–30% discrepancy applies only to a narrow band near the boundary

**Required revision:**
- Qualify: *"For the $\sim$25% of points that disagree by more than 1%, the
  discrepancy reaches 20–30% near the metastability boundary"*
- Do not use "catastrophically" without quantification
- Add a fractional-error map (see MISSING_ANALYSIS.md M4)

---

## Impact 6: "c₆ rescues $\sim$110 points corresponding to 0.3–0.5 GeV in $M_t$" (Abstract, §V)

**Severity: MAJOR** (would be FATAL if not already flagged in MISSING_ANALYSIS.md)

**Where it appears:**
- Abstract lines 46–49: "rescues a continuous strip... approximately 110 grid
  points corresponding to $\sim$0.3–0.5~GeV in $M_t$"
- §V.B lines 755–761 and Fig. 7 caption

**What the code/data show:**
- The c₆ sweep solver (`archive/solver_c6_sweep.cpp`) is a **separate executable**
  not wired into the pipeline
- **No c₆-column data exists** in `results/`
- No run script for c₆ sweep exists
- The reported numbers (110, 0.3–0.5) are hardcoded in the LaTeX with no
  supporting data file

**Required revision:**
- If c₆ sweep was actually performed on a different system, the data must be
  added to the repository and the analysis reproducible
- If not, these claims **must be removed** from the abstract and main text
- At minimum, the paper must state: *"The c₆ sweep was performed on a separate
  compute cluster; the reproducibility data is available from the authors upon
  request"*
- The abstract should not feature this as a headline result without supporting
  data

---

## Impact 7: "250,000 evaluation points" (Abstract, §III.E)

**Severity: MINOR**

**Where it appears:**
- Abstract line 29: "computed across 250,000 parameter points"
- §III.E line 342: "a total of 250,000 evaluation points"

**What the code/data show:**
- The CSV contains **1,737,982 points**
- The 250k refers to a subset of the grid (500×500) in the range
  M_t ∈ [155, 185], M_h ∈ [115, 135]
- The full scan extends to M_t ∈ [0, 250], M_h ∈ [0, 250] at 0.25 GeV spacing

**Required revision:**
- Clarify: *"The primary analysis uses a 500×500 uniform grid (250,000 points)
  covering the physically relevant region, embedded in a larger 1001×1001 grid
  (1,000,000+ points) for full coverage"*
- This is a minor inaccuracy that undersells the computation

---

## Impact 8: "Fubini-Lipatov breakdown creates a sandwiched unstable region" (§IV.D)

**Severity: FATAL**

**Where it appears:**
- §IV.D header and all text in that subsection
- Abstract line 52–53 (implicitly, through the claim of numerical necessity)

**What the code/data show:**
- Please see CLAIM_AUDIT.md Claim 3 for the data analysis
- No sandwiched (2→3→2) pattern exists in `results/numerical_data.csv`
- Negative S_exact occurs predominantly at high M_t > 231, not at low M_h
- S_kin and S_pot are not exported, so the virial ratio cannot be verified

**Required revision:**
- Either remove §IV.D entirely, or
- Substantially rewrite it to report: (a) the absence of a sandwiched pattern
  in the data, (b) the distribution of negative S points, and (c) that the
  Fubini ansatz breakdown is the explanation for negative S, not a physical
  sandwiched region
- The term "sandwiched" must be replaced with "Fubini-Lipatov ansatz failure regime"

---

## Impact 9: "The true bounce in this regime... yields a genuine, small positive action" (§IV.D, line 751)

**Severity: MAJOR**

**Where it appears:**
- §IV.D line 751: "the true bounce in this regime — obtainable by direct
  numerical integration of the O(4)-symmetric shooting equations with the full
  running potential — yields a genuine, small positive action"

**What the code shows:**
**This computation was never performed.** The paper asserts a result it did not
compute. The statement is speculative.

**Required revision:**
- Remove this sentence, or
- Replace with a plans-verb: *"A dedicated ODE bounce solver would be required
  to determine the true action in this regime; such a computation is beyond the
  scope of this paper"*

---

## Summary

| # | Claim | Severity | Requires |
|---|-------|----------|----------|
| 1 | "Exact numerical bounce solver" | **FATAL** | Language replacement throughout |
| 2 | "Numerical errors below $10^{-11}$" | **FATAL** | Add ansatz systematics to error budget |
| 3 | "Adaptive RK4 step sizing" | **MAJOR** | Correct to fixed-step description |
| 4 | "EOM residuals validate profiles" | **MAJOR** | Reframe as ansatz accuracy diagnostic |
| 5 | "20–30% error near boundary" | **MAJOR** | Qualify with quantified fractions |
| 6 | "c₆ rescues 110 pts" | **MAJOR** | Provide data or remove |
| 7 | "250,000 evaluation points" | MINOR | Clarify |
| 8 | "Sandwiched unstable region" | **FATAL** | Remove or rewrite entirely |
| 9 | "True bounce gives small positive action" | **MAJOR** | Remove uncomputed speculation |

**Minimum revision required before resubmission:**
1. Rewrite all instances of "exact numerical solver" (Impacts 1, 2, 3)
2. Remove or fully rewrite the "sandwiched region" section (Impact 8)
3. Provide c₆ data or remove the claim (Impact 6)
4. Remove uncomputed speculation about the true bounce (Impact 9)
5. Add ansatz systematic error to Table I (Impact 2)

**Estimated revision depth:** The abstract, introduction, methodology (§III),
results (§IV.D), and error budget (§IV) all require significant rewriting.
Approximately 40–60% of the paper text must change.
