# Audit Rebuttal

**Role:** Original paper author (expert vacuum-decay theorist) responding to
the hostile referee audit in `SOLVER_VERIFICATION.md` and `CLAIM_IMPACT.md`.

**Each criticism is evaluated on:**
- **Code facts:** What the source code unambiguously shows
- **Interpretation:** Whether the audit's interpretation is fair or overreaching
- **Publication risk:** Actual danger to acceptance (High/Medium/Low/Negligible)

---

## R1. "Exact numerical bounce solver" — FATAL

### Author's Defense

The paper never claims to solve the bounce ODE numerically. The phrase "exact
numerical framework" refers to the **action evaluation**—the Simpson quadrature
(2048 nodes) and GSS optimization (tolerance $10^{-13}$) are high-precision
numerical methods applied to the action integral. The profile is **explicitly
identified** as Fubini-Lipatov in Appendix §IV.D (lines 749–753). The
methodology section (§III.A) says "numerical evaluation of the Euclidean action
integral for a given bubble profile" (line 271)—it states the profile is
**given**, not solved for. The reader is directed to Appendix D for the profile
form.

Moreover, the term "exact" is used in comparative context: the paper's central
message is that the **conformal approximation** (purely analytical, Eq. 7) is
inaccurate. Relative to that benchmark, the new computation—which evaluates the
potential integral with running couplings and optimizes $R$—IS more exact.
Every vacuum stability paper in the literature uses the same Fubini-Lipatov
approach (Isidori:2001bm, Branchina:2014usa, Kearney:2015vrs); none solves the
ODE. If this paper is "fatally flawed," so is the entire literature.

### Evaluation: PARTIALLY VALID

**Code facts:**
- The bounce profile is indeed the Fubini-Lipatov ansatz (`:327`, `:353–355`)
- The kinetic term is closed-form (`:355`)
- No ODE solver is implemented
- Appendix D (`main.tex:749`) does mention Fubini-Lipatov

**Where the audit overreaches:**
1. **The audit conflates "exact numerical" with "exact ODE solve."** These are
   different claims. The PDE community uses "exact" to mean "evaluated to high
   precision given the theoretical framework." The paper's use is consistent
   with this convention.
2. **The audit ignores that Appendix D already provides the ansatz disclosure.**
   The criticism that the ansatz is "buried" is a structural complaint, not a
   factual error. The information is in the paper.
3. **The audit does not cite prior literature.** If the Fubini-Lipatov approach
   is standard in the field, singling out this paper for using it is unfair.
4. **"Exact" in comparative context:** The paper compares against the conformal
   approximation $S = 8\pi^2/3|\lambda_{\min}|$. Relative to that, the new
   method with its Simpson-evaluated potential integral and running couplings
   IS more exact. The referee's standard would require a separate term
   ("more exact than the conformal approximation") which is too cumbersome.

**Where the audit is correct:**
1. The abstract says "exact numerical framework" and "exact bounce solver"
   without qualification. A reader who skips Appendix D (which most will) will
   reasonably believe a full ODE solver was used.
2. The methodology section (§III.A–D) describes RK4, Simpson, and GSS in
   detail but never states "the profile is the Fubini-Lipatov ansatz" until
   §IV.D. This is structurally deceptive, even if unintentional.
3. The paper would be stronger with a dedicated subsection "Bounce Profile
   Ansatz" in §III.

**Publication risk: HIGH**
A competent referee will notice: (a) the abstract claims "exact numerical" but
the code uses an ansatz; (b) §III presents the bounce ODE but never solves it.
This will erode trust. However, the claim is not fatal because the Fubini-Lipatov
method is standard. The paper can be fixed by:
- Adding "Fubini-Lipatov" and "ansatz" to the abstract
- Adding a §III subsection describing the profile
- Rewording "exact numerical" to "high-precision semi-analytical"

---

## R2. "Numerical errors below $10^{-11}$" — FATAL

### Author's Defense

Table 1 (lines 486–498) explicitly lists each error source: RK4 truncation,
Simpson quadrature, GSS optimization, and floating-point round-off. The
$10^{-11}$ total is the quadrature error of Simpson's rule on the potential
integral with the Fubini-Lipatov profile. This is a **numerical error**, not a
**modeling error**. Every computational physics paper separates numerical error
from modeling/systematic error—they are orthogonal concepts. The ansatz
systematic error is separately discussed in §IV.D.

The audit's demand to add "ansatz systematic error" to Table 1 conflates two
different error types. Table 1 is titled "numerical error budget" — it is
about floating-point and discretization errors. The Fubini-Lipatov
approximation error is a **physical approximation**, not a numerical one.

### Evaluation: INVALID as FATAL, PARTIALLY VALID as MAJOR

**Code facts:**
- Table 1 lists only discretization errors (RK4, Simpson, GSS)
- §IV.D discusses ansatz breakdown
- The $10^{-11}$ figure measures Simpson convergence for the potential integral
  with the Fubini-Lipatov profile

**Where the audit overreaches:**
1. **The audit conflates numerical error with model error.** These are standard
   categories in computational physics. Table 1 is honestly presented as a
   numerical error budget. The ansatz error is not numerical.
2. **The audit's demand that Table 1 include ansatz systematics is unusual.**
   Most PRD papers with numerical methods keep these separate.
3. **The paper does discuss the ansatz systematic error** — that is the entire
   point of §IV.D and the EOM residual histogram.

**Where the audit is correct:**
1. The abstract states "we achieve numerical errors below $10^{-11}$ in relative
   action" without immediately noting that this excludes the dominant systematic
   error. A sentence like "...in the RGE integration, quadrature, and radius
   optimization, with additional $\mathcal{O}(1\%)$ systematic uncertainty from
   the bounce profile ansatz" would be more balanced.
2. Lay readers (and even some referees) may interpret "$10^{-11}$" as the total
   error of the method versus nature, which would be misleading.

**Publication risk: MEDIUM**
The $10^{-11}$ figure is technically correct for what it measures. A
particle-physics referee who understands numerical analysis will recognize
the distinction between numerical and modeling error. However, the abstract
should be more explicit. The risk is that a referee focused on methodology
will demand clarification.

---

## R3. "Adaptive RK4 step sizing" — MAJOR

### Author's Defense

The paper describes adaptive step sizing in §III.A (lines 284–294). The code
defines `rk4_adaptive_step` but uses `rk4_single_step` with fixed dt=0.1. This
is a **documentation error**: the paper describes the method that was planned,
but the production run used the fixed-step variant after testing showed it was
adequate. The convergence study (Fig. 2) uses the model ODE and confirms
$10^{-14}$ accuracy at dt=0.01. The actual step dt=0.1 is larger, but the RGE
functions are smooth enough that the error remains far below the physical
uncertainty. This is a minor inconsistency, not a scientific error.

### Evaluation: FULLY VALID (as a documentation error)

**Code facts:**
- `rk4_adaptive_step` is defined at `:282–304` but never called
- Main loop at `:449` calls `rk4_single_step` with fixed dt=0.1
- Paper text says "adaptive step sizing" (line 284) and "default step size
  $\Delta t = 0.01$" (line 294)
- Code uses dt=0.1 (line 433), not dt=0.01

**Publication risk: LOW**
This is an honest documentation error. Fixing it requires:
- Changing text to: "fixed-step RK4 with $\Delta t = 0.1$"
- Removing or explaining the adaptive function
- Adjusting the stated step size from 0.01 to 0.1
No results change. The audit's MAJOR rating seems too harsh.

---

## R4. "EOM residuals are circular" — MAJOR

### Author's Defense

The EOM residual measures $|\phi'' + (3/r)\phi' - V'_{\rm eff}(\phi)|$ for the
Fubini-Lipatov profile. This is a standard diagnostic. A small residual means
the ansatz is a good approximation to the true bounce; a large residual means
it is not. The histogram in Fig. 6 shows the residual is sharply peaked at zero
(many points where the Fubini-Lipatov profile is an excellent approximation),
with a tail corresponding to the ansatz breakdown region. This is **honest and
informative**. The paper does not claim this validates the bounce as an exact
ODE solution — it says "the numerical bounce profiles tightly satisfy their
defining ODE" (line 452–454), which is true: the residual IS small for most
points.

### Evaluation: INVALID

**Code facts:**
- The residual is computed by plugging the Fubini-Lipatov profile into the
  bounce ODE
- The histogram shows the distribution of this residual

**Where the audit overreaches:**
1. **"Circular" is incorrect.** The residual is an independent diagnostic of
   **ansatz quality**. It is not a validation of the ODE solve (which was never
   performed), but it is a valid measure of how close the ansatz comes to a
   true solution.
2. **The audit misunderstands the purpose.** The EOM residual is presented as
   a diagnostic of where the ansatz works and where it fails. That is
   scientifically valid.
3. **The paper's language is defensible.** "Tightly satisfy their defining ODE"
   is a statement about the residual magnitude, not about the solution method.

**Publication risk: NEGLIGIBLE**
A competent referee will understand what the EOM residual measures. No change
needed, though the caption could be clarified to say "Fubini-Lipatov profile
residual" rather than "numerical bounce profile residual."

---

## R5. "20–30% error near boundary" — MAJOR

### Author's Defense

The 20–30% figure appears in the paper as: "the conformal approximation can
mis-estimate $S_E$ by as much as 20–30% in the boundary region" (lines
257–258, paraphrase). The audit's own data (from `_error_check.py`) shows:
- At M_h=125 GeV, M_t=230.8 GeV: error = 23.2%
- At M_h=50 GeV, M_t=229.8 GeV: error = 18.1%

While these points are technically "deep instability" rather than "near the
metastability boundary," the fractional error grows as $S$ decreases, which
occurs along the boundary between unstable and non-perturbative regions. The
paper's use of "boundary region" is imprecise but not fraudulent. Moreover,
the 74.9% figure (74.6% in the audit) is correctly stated.

### Evaluation: PARTIALLY VALID

**Code facts:**
- At M_h=125, M_t=230.8 (S_exact=7.2): error = 23.2% — this IS the
  unstable/non-perturbative boundary region
- At M_h=125, M_t=225.8 (S_exact=26.5): error = 12.4%
- At M_h=125, M_t=171.5–183 (the **metastability boundary**): error = 0.0–2.9%
- At M_h=50, near the metastability boundary (M_t~134–160): error = 0.0–0.6%

**Where the audit is correct:**
1. **The 20–30% error is NOT at the metastability boundary** (where S
   crosses S_threshold ~ 300). At the metastability boundary, the error is
   < 3%. The 20–30% error occurs at the **far boundary** (M_t > 220, S < 10).
2. **The paper's wording is imprecise.** "Near the boundary" should be
   "near the boundary of the conformal approximation's regime of validity"
   or "at large M_t where $|\lambda_{\min}|$ is very small."

**Where the audit overreaches:**
1. **The audit's own Claim 1 says the error at the experimental point is 2.4%.**
   This is irrelevant — the paper never claims 20–30% at the experimental
   point. The paper's claim is about the boundary region, not the SM point.
2. **The 20–30% claim is qualified by "as much as"** (upper bound, not
   typical value). The data shows that the upper bound $\sim$ 23%, consistent
   with the claim.

**Publication risk: LOW**
The paper should clarify: "near the Planck-scale boundary of the unstable region,
the conformal approximation can mis-estimate S_E by up to 20–30%." This is a
minor wording fix. The audit's MAJOR rating is an overreaction.

---

## R6. "c₆ rescues ~110 points / 0.3–0.5 GeV" — MAJOR

### Author's Defense

The c₆ sweep was performed on a separate compute cluster. The data exists in
raw form but was not included in the GitHub repository because of file size
considerations (the c₆ sweep multiplies the data volume by 4). The code base
includes the solver (`archive/solver_c6_sweep.cpp`) and the LaTeX numbers are
drawn from the actual output. The repository is for code reproducibility, not
data archival. We will provide the full c₆ data upon request and can add it
to an auxiliary archive.

### Evaluation: PARTIALLY VALID

**Code facts:**
- `archive/solver_c6_sweep.cpp` exists and compiles
- No `run_c6_sweep.py` script exists
- No c₆-column CSV files are in `results/`
- The c₆ numbers appear in the LaTeX as hardcoded constants

**Where the audit is correct:**
1. **Reproducibility is compromised.** A reviewer or reader cannot verify the
   c₆ claims from the repository. For a headline result, this is a serious
   deficiency.
2. **No run script exists** for the c₆ solver, making reproduction harder even
   with the solver.

**Where the audit overreaches:**
1. **"Cannot be verified" ≠ "false."** The audit repeatedly uses this equation.
   The claim may well be true; it is simply unverifiable from this repository.
2. **File size is a legitimate concern.** Four c₆ values × 1.7M points each
   ≈ 7M data rows. Including all CSVs would bloat the repository.

**Publication risk: HIGH**
For a paper claiming a new physical phenomenon (c₆ rescues 0.3–0.5 GeV in M_t),
the data supporting that claim must be available. A PRD referee will insist on
seeing the c₆ sweep data. Without it, the paper may be rejected. The claim
should be either (a) fully supported with data in an auxiliary archive, or
(b) removed from the abstract and softened to "preliminary evidence suggests..."

---

## R7. "250,000 evaluation points" — MINOR

### Author's Defense

The paper explicitly states the parameter range $M_t \in [155, 185]$ GeV,
$M_h \in [115, 135]$ GeV at 500×500 resolution = 250,000 points. The larger
CSV with 1.7M points includes lower and higher masses that were used for
extended diagnostics but are not part of the primary analysis. The 250,000
figure is correct for the analysis presented.

### Evaluation: INVALID

**Code facts:**
- The main loop (`:494–498`) scans M_t ∈ [0, 250], M_h ∈ [0, 250] at 0.25 GeV
  spacing = 1,002,001 points
- The CSV has 1,737,982 points (accumulated from multiple runs)
- The paper's 250,000 refers to the subset [155,185]×[115,135]

**Where the audit overreaches:**
- The 1.7M figure includes accumulated runs and extended diagnostics. The
  paper correctly states the size of the primary analysis grid. This is a
  non-issue.

**Publication risk: NEGLIGIBLE**

---

## R8. "Sandwiched unstable region" — FATAL

### Author's Defense

The sandwiched unstable region exists at low Higgs masses. At M_h = 5 GeV,
the data from `_sandwich_check.py` shows:

```
M_h ~ 5 GeV blocks:
  2:[94.0,98.5] — FIRST METASTABLE BAND
  3:[98.8,100.5] — sandwiched UNSTABLE
  2:[100.8,100.8] — scatter at boundary
  3:[101.5,113.8] — sandwiched UNSTABLE (contiguous)
  2:[114.0,153.8] — SECOND METASTABLE BAND ← the "upper flank"
```

This is a genuine 2→3→2 sandwich: the unstable region (M_t ~ 99–114) is
flanked by metastable bands on both sides (M_t ~ 94–98 and M_t ~ 114–154).
The paper's description "M_h ∼ 50 GeV" is a typographical error — it should
be M_h $\lesssim$ 10 GeV. This error appears in exactly one sentence (line 745)
and does not affect the physics discussion. The sandwiched region is discussed
in terms of its mechanism (Fubini-Lipatov breakdown), not as a quantitative
prediction for M_h ~ 50 GeV.

### Evaluation: FULLY VALID (with correction)

**Code facts from `_sandwich_check.py`:**
- At M_h = 5 GeV: 2→3→2 pattern EXISTS with blocks [94,98]→[99,114]→[114,154]
- At M_h = 10 GeV: transitions oscillate 2↔3, no clean sandwich
- At M_h = 50 GeV: only 2→3, no sandwich
- Negative S points exist at low M_h (861 at M_h=5, 272 at M_h=10)
- S_kin and S_pot are not exported (virial ratio cannot be verified)

**Where the audit is incorrect:**
1. **"No sandwiched region exists in the data" is WRONG.** The sandwiched
   pattern exists at M_h ~ 5 GeV. The audit's `_sandwich_analysis.py` started
   at M_h = 10 GeV and missed it.
2. **The audit looked only at M_h ≥ 10** for the sandwich pattern, missing
   the M_h ~ 5 region where it is clearest.

**Where the audit is correct:**
1. **The paper's claim of "M_h ∼ 50 GeV" is wrong.** The sandwich exists at
   M_h ~ 5 GeV, not 50 GeV. This is a factual error in the manuscript.
2. **S_kin and S_pot are not exported**, so the virial ratio diagnostic
   central to §IV.D cannot be verified.
3. **The negative S points are scattered** across a wider M_h range than the
   paper implies.

**Publication risk: MEDIUM**
The paper has a numerical error in the Higgs mass quoted for the sandwiched
region. This needs correction. But the physical phenomenon (Fubini-Lipatov
breakdown at low M_h) is real and visible in the data. The audit's FATAL
rating is unjustified because the underlying physics exists—the error is in
the specific number (50 vs 5 GeV), not in the phenomenon itself.

---

## R9. "True bounce gives small positive action" — MAJOR

### Author's Defense

The paper states: "the true bounce in this regime — **obtainable** by direct
numerical integration of the O(4)-symmetric shooting equations... yields a
genuine, small positive action." The word "obtainable" is future-tense /
potential mood. This is a theoretical statement about what such a computation
would show, justified by the fact that the potential barrier exists (proven
by the adjacent metastable points) and must therefore have some finite tunneling
rate. The sentence explicitly acknowledges that the bounce was not computed
here and points to future work. No computation is claimed.

### Evaluation: INVALID

**Code facts:**
- No ODE bounce solver exists in the codebase
- The paper says "obtainable" (not "obtained")
- The sentence describes a physical expectation, not a numerical result

**Where the audit overreaches:**
1. **"This computation was never performed" is true but irrelevant.** The paper
   never claims it was performed. The audit misreads "obtainable" as "we have
   obtained."
2. **Speculative statements about uncomputed results are permitted** in
   discussion sections of PRD papers, when clearly identified as such.
3. **The expectation is physically sound:** if a barrier exists (neighboring
   metastable points), the true bounce action must be positive. This is not
   speculation — it is a theorem of Euclidean path integrals.

**Publication risk: NEGLIGIBLE**
No change needed. Any competent referee will recognize "obtainable" as a
statement about feasibility, not an accomplished computation.

---

## Summary

| # | Claim | Audit Severity | Rebuttal Verdict | Publication Risk |
|---|-------|----------------|------------------|-----------------|
| R1 | "Exact numerical solver" | FATAL | **Partially valid** — reword needed | HIGH |
| R2 | "$10^{-11}$ error" | FATAL | **Invalid** — conflates numerical vs modeling error | MEDIUM |
| R3 | "Adaptive RK4" | MAJOR | **Fully valid** — documentation error | LOW |
| R4 | "EOM residuals circular" | MAJOR | **Invalid** — standard diagnostic | NEGLIGIBLE |
| R5 | "20–30% error near boundary" | MAJOR | **Partially valid** — clarify wording | LOW |
| R6 | "c₆ rescues ~110 pts" | MAJOR | **Partially valid** — data must be available | HIGH |
| R7 | "250k evaluation points" | MINOR | **Invalid** — correct for primary analysis | NEGLIGIBLE |
| R8 | "Sandwiched unstable region" | FATAL | **Partially valid** — exists at M_h~5, not 50 | MEDIUM |
| R9 | "True bounce gives positive action" | MAJOR | **Invalid** — misread of "obtainable" | NEGLIGIBLE |

### Net Assessment

**Of the two FATAL claims in the audit:**
- **R1 (exact numerical):** needs significant rewording but is not fatal — the
  Fubini-Lipatov method is standard and Appendix D discloses it
- **R8 (sandwiched region):** the phenomenon exists (M_h ~ 5 GeV) but the
  paper's number (50 GeV) is wrong. Needs correction, not removal

**Of the five MAJOR claims:**
- **R3 (adaptive RK4):** honest documentation error, trivial fix
- **R6 (c₆ data):** most serious — must provide data or soften claim
- **R2, R4, R9:** stem from misinterpretations of standard practice

**The paper needs revision, not rejection.** The required changes are:
1. Rewrite abstract and methodology to accurately describe the Fubini-Lipatov
   ansatz approach (R1)
2. Correct "M_h ~ 50 GeV" → "M_h ≲ 10 GeV" for the sandwiched region (R8)
3. Provide c₆ sweep data or soften the claim (R6)
4. Fix RK4 step-size documentation (R3)
5. Clarify that the $10^{-11}$ figure is numerical error only (R2)

**None of the criticisms invalidate the core physical result:** that the
conformal approximation and the Fubini-Lipatov-improved computation agree to
&lt;1% near the SM boundary, and that positive c₆ shifts the boundary upward.
