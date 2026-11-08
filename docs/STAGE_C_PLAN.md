# Stage C Plan — Architecture and Numerical-Method Comparison

**Project:** SMVac — Naman Goyal (namangoyal31), IISc Bangalore
**Status:** plan agreed with the owner after the Stage A/B review (`f3d2e25`). This document
records the decisions and the proposed method comparison *before* implementation. The evidence
for the final choices is collected in `docs/NUMERICAL_METHOD_COMPARISON.md` and
`docs/VALIDATION_REPORT.md`.

---

## 1. Central scientific question

> How accurately do commonly used analytical approximations describe the underlying numerical
> problem, where do they fail, and how much can we improve upon them?

The project is organised as:

1. reproduce the literature approximation;
2. state its assumptions and domain of validity;
3. build and validate the most reliable practical numerical reference;
4. quantify the discrepancy;
5. test improved analytical estimates and bounds;
6. demonstrate everything with convergence studies, error analysis and reproducible figures.

## 2. The stated problem

**Problem P (tree-level semiclassical decay of the RG-improved potential).** For SM inputs
$(M_h,M_t,\alpha_3(M_Z))$:

1. Compute the MS-bar couplings at $\bar\mu=M_t$ (Buttazzo et al. 2013, eqs. 55–60) and evolve
   them with the 3-loop SM RGEs (their App. B) to $M_{\rm Pl}=1.22\times10^{19}$ GeV.
2. Build the RG-improved quartic $\lambda_P(\phi)$ at $\bar\mu=\phi$ for a prescription $P$:
   * `tree`: the MS-bar $\lambda(\bar\mu=\phi)$;
   * `cw1`: the one-loop Coleman–Weinberg $\lambda_{\rm eff}$ of §4 of the methodology, with the
     corrected Z coupling and no $e^{4\Gamma}$.

   Couplings are frozen outside $[\max(M_t,v),M_{\rm Pl}]$.
3. Absolute stability holds iff $\min_\phi\lambda_P\ge0$.
4. Otherwise compute the O(4) bounce action $S$ of $V=\lambda_P(\phi)\phi^4/4$ in flat space,
   its scale $\Lambda_B=1/R$, and the past-light-cone decay probability
   $\wp=0.15\,\Lambda_B^4e^{-S}/H_0^4$ (Buttazzo et al. eqs. 80–82). The vacuum is metastable
   iff $\wp<1$.

The φ⁶ term is **not** part of P. It exists only as an explicit, opt-in deformation
`bsm.phi6` $=c_6\phi^6/(6\Lambda^2)$ that is recorded in every output's metadata and never used
in SM results or figures.

**Not included** (stated as limitations): the one-loop fluctuation determinant, gravitational
corrections, the two-loop $\lambda_{\rm eff}$ and $e^{4\Gamma}$, and thermal effects.

## 3. Vocabulary (used consistently in code and documents)

| Term | Meaning | Examples |
|---|---|---|
| **numerical reference** | a numerical solution of P with a quantified numerical uncertainty $\delta S$; never called "exact" | the chosen bounce solver |
| **analytical approximation** | closed-form or quadrature formula with no error guarantee | $S_{\rm FL}=8\pi^2/(3\lvert\lambda_{\min}\rvert)$ (literature), the curvature-corrected $S_{\rm curv}$ (§5) |
| **rigorous bound** | an inequality proven under stated assumptions | $S\ge S_{\rm FL}$ (Pohozaev + sharp Sobolev) |
| **variational estimate** | an upper estimate that is a rigorous bound only under an unproven hypothesis | $S_{\rm CGM}=8\pi^2/(3\lvert\min_h\bar\lambda(h)\rvert)$ (needs CGM attainment) |
| **legacy** | original implementation, kept for traceability | `S_exact` |

## 4. Domain of validity (owner decision 4)

* **Controlled region** (headline figures and conclusions): $M_h\in[120,130]$ GeV and
  $M_t\in[166,182]$ GeV, with $\alpha_3(M_Z)=0.1184$ by default.
  * The matching formulas are first-order expansions about (125.15, 173.34) GeV. Quadratic
    terms, estimated from tree-level relations, stay below the quoted theory uncertainty
    ($\sim3\times10^{-4}$ on $\lambda$) inside this window.
  * Points whose bounce probes field values within a decade of $M_{\rm Pl}$, or whose
    $\lambda_{\min}$ sits on the frozen Planck boundary, are flagged **Planck-sensitive**
    and shown hatched.
* **Exploratory region** $[0,250]^2$ GeV: recomputed with the new code at coarse resolution
  and labelled "uncontrolled extrapolation" on every figure and in every output header. The
  historical legacy outputs are preserved in `legacy/data/`.

## 5. Analytical approximations and bounds to compare

| Name | Formula | Status |
|---|---|---|
| $S_{\rm FL}$ (literature) | $8\pi^2/(3\lvert\lambda_{\min}\rvert)$ | rigorous **lower bound** for every finite-action bounce |
| $S_{\rm curv}$ (new) | $8\pi^2\big/\big[3(\lvert\lambda_{\min}\rvert-\tfrac{13}{72}\lambda''_{\min})\big]$, with $\lambda''=d^2\lambda/d(\ln\phi)^2$ | analytical approximation: first-order perturbation theory about the FL instanton at the minimum with a quadratic $\lambda$; the 13/72 is $\tfrac12\mathrm{Var}(\ln w)$ under the FL weight $6w(1-w)$ |
| $S_{\rm CGM}$ (new) | $8\pi^2/(3\lvert\min_h\bar\lambda(h)\rvert)$, $\bar\lambda(h)=6\int_0^1w(1-w)\lambda(hw)dw$ | variational upper estimate |
| `S_exact` (legacy) | FL first-order estimate with $\lambda_c=\lambda(1/R)$ and φ⁶ | legacy baseline, known unsound |

**Expected behaviour (hypothesis to test).**

* Deep in the metastable region, $S_{\rm FL}$ errs by $\mathcal O(10^{-3})$, and $S_{\rm curv}$
  and $S_{\rm CGM}$ reduce this to $\mathcal O(10^{-6})$.
* Near the absolute-stability boundary, the potential develops a second, nearly degenerate
  minimum below $M_{\rm Pl}$. The bounce becomes thin-walled and $S$ grows much faster than
  $1/\lvert\lambda_{\min}\rvert$. We therefore expect $S_{\rm FL}$, and the FL-shaped trial
  inside $S_{\rm CGM}$, to fail badly there.

This is the regime that matters for "near-criticality".

## 6. Proposed numerical-method comparison (owner decision 1)

All methods solve the same problem P and use the same $\lambda_P$ representation, so that only
the method differs. The variables used below are $\tau=\ln r$ and $f=r\phi$, giving
$f''=f+f^3(\lambda+\tfrac14\beta)$ with $\beta=d\lambda/d\ln\phi$ (METHODOLOGY §6.1).

### M0. Legacy FL-ansatz evaluation with golden-section search (baseline)

* **Why included:** traceability.
* **Expected:** $\mathcal O(1)$ errors, negative values, sentinels. It is not a solution of the
  bounce equation at all.

### M1. Classical overshoot/undershoot shooting in physical variables $(r,\phi)$

Rescaled units, a series start at $r\to0$, adaptive Dormand–Prince integration, and bisection
on $\phi(0)$. This is the textbook method, so it serves as the community-standard numerical
comparator.

* **Expected weaknesses:**
  * the $3\phi'/r$ singularity at the origin;
  * $\phi=0$ is a *non-hyperbolic* (massless) end point, approached algebraically
    ($\phi\sim C/r^2$), which makes over/undershoot detection slow and truncation errors
    algebraic;
  * a field range of up to $10^{19}$ GeV must be handled by rescaling.

### M2. Shooting in scale-covariant variables $(\tau,f)$

The same bisection logic as M1, but $f=0$ is a *hyperbolic* saddle with known exponents
$e^{\pm\tau}$.

* **Why it should be more accurate:**
  * there is no origin singularity: the start uses the exact linearised solution;
  * $f=O(|\lambda|^{-1/2})$ has an $O(1)$ dynamic range;
  * the energy $E=\tfrac12f'^2-\tfrac12f^2-\tfrac14\lambda f^4$ obeys $dE/d\tau=\tfrac14\beta f^4$,
    giving a sharp over/undershoot criterion and two built-in identities (Derrick, dilatation);
  * the action is stationary, so it inherits only the *square* of the shooting error.
* **Expected weakness:** thin walls. When the bounce sits near a true vacuum for many e-folds,
  the required precision of $\phi(0)$ grows exponentially. Bisection in double precision then
  fails beyond roughly 35 e-folds.

### M3. Global collocation BVP (4th-order Lobatto IIIA, adaptive mesh; `scipy.solve_bvp`) on a truncated $\tau$ interval

Robin conditions $f'=\pm f$ at the ends (exact asymptotics), with an FL initial guess.

* **Why it may be more reliable:** a global Newton solve has no exponential sensitivity, so it
  should handle thin walls that defeat shooting.
* **Expected weakness:** the dilatation quasi-zero mode makes the Jacobian nearly singular when
  the running is slow, so convergence depends on the initial guess. Accuracy is limited by the
  4th-order mesh.

### M4. Chebyshev spectral collocation with Newton iteration

Same truncated domain and Robin conditions as M3, with a smooth Chebyshev representation of
$\lambda(\ln\phi)$.

* **Why it may be the most accurate:** exponential convergence in the number of nodes for an
  analytic solution. It is independent of the shooting machinery, so it makes a genuine
  cross-check.
* **Expected weakness:** the same soft mode as M3. Note, though, that the *action* is stationary
  along that mode, so its error is second order.

### M5 (optional). Constrained minimisation (Coleman–Glaser–Martin reduced problem: minimise $T$ at fixed $V_{\rm int}$)

* **Why:** it turns the saddle-point problem into a minimisation (robust; the idea behind
  gradient-flow solvers).
* **Expected weakness:** the soft mode slows gradient methods. It will be tried only if M3/M4
  prove fragile in the thin-wall regime.

### Benchmarks

* **B1:** constant $\lambda$, where $S=8\pi^2/(3|\lambda|)$ exactly.
* **B2:** a parabolic $\lambda(x)$ in the thick-wall regime (smooth, no closed form).
* **B3:** a model with a second minimum and tunable degeneracy (thick → thin wall).
* **B4:** SM points in the controlled region: the measured point, deep metastable, near the
  metastability boundary, and a sequence approaching the stability boundary.
* **B5:** exploratory points (low $M_h$, large $M_t$).

### Metrics

* accuracy against B1;
* self-convergence under refinement (tolerance, mesh/N, potential-table resolution);
* Derrick and dilatation residuals;
* cross-method agreement;
* wall-clock cost per point;
* robustness: failure rate over a random sample of the domain, especially near the stability
  boundary;
* sensitivity to initial guesses and tolerances.

### Selection rule

The primary reference is the method with the smallest *demonstrated* uncertainty at acceptable
cost and with the best robustness. A second, structurally different method is kept as the
cross-check. A hybrid (shooting for thick walls, collocation for thin walls) is acceptable if
the evidence favours it.

**Uncertainty assigned to the reference:** $\delta S$ = max(refinement change, cross-method
difference, identity residual × S, change under a 2× change of the potential-table resolution).

## 7. Software architecture (owner decision 5)

```
CMakeLists.txt            C++14 library + apps + CTest
include/smvac/  src/      constants, couplings, beta_functions (Buttazzo App. B), matching
                          (buttazzo2013 | legacy), rge (RK4 + Hermite dense output),
                          potential (tree | cw1; opt-in bsm.phi6), analytic (FL, curv, CGM),
                          bounce (chosen numerical reference + diagnostics + uncertainty),
                          decay (eqs. 80-82), stability (classification with explicit statuses),
                          numerics/ (RK4, Dormand-Prince, Brent, adaptive quadrature, Hermite)
apps/                     smvac_point (full JSON report), smvac_scan (grid -> CSV + metadata)
python/smvac/             config, parallel orchestration, I/O, plotting
python/smvac/reference/   independent Python implementations (DOP853 RGE, M1-M4) for validation
configs/                  scan definitions (controlled, exploratory, benchmarks)
validation/               method comparison, convergence studies, C++-vs-Python cross-checks
tests/                    unit, physics, convergence and regression tests (independent references only)
results/  figures/        curated outputs with metadata; regenerated figures
legacy/                   original code, scripts, data, figures (git mv; frozen; README with known defects)
audit/                    Stage A/B tooling
docs/                     audit, methodology, reviews, method comparison, validation report, reproducibility;
                          docs/archive/ for predecessor planning documents
```

* **Language split.** C++ for the production scans (thousands of bounce solves), Python for the
  independent reference implementations and plotting. Agreement between two implementations in
  two languages is itself a validation result.
* **Compiler constraint:** C++14, because the available MinGW.org GCC 6.3 lacks C++17 library
  features. Parallelism is by process (orchestrated from Python), as before.

## 8. Physics corrections in the new core

Each correction lands in a separate, validated commit:

1. matching to Buttazzo eqs. 55–60 (configurable inputs; legacy set kept as an option);
2. the three RGE transcription errors;
3. the Z coupling in $\lambda_{\rm eff}$;
4. φ⁶ removed from P (opt-in deformation only);
5. the literature decay criterion with $\Lambda_B=1/R$ from the bounce;
6. Brent root-finding and minimisation for the instability scale and $\lambda_{\min}$;
7. explicit statuses for stable, non-perturbative, Planck-sensitive, no-bounce and solver
   failure, so that failures are never silent.

## 9. Repository clean-up (owner decision 5, with safeguards)

* `git mv` the monolith, `scripts/`, `SMVacuumDecay/`, `data/` and the old `figures/` into
  `legacy/` (no history rewriting, no extra storage), with `legacy/README.md` documenting
  provenance and defects.
* Audit each CSV before moving.
  * `data/numerical_data.csv`: historical final output of the legacy wide scan (with
    duplicates, a corrupted row, and reproducible unique content).
  * `data/analytical_data.csv`: historical output of `classify_buttazzo`.
  * `SMVacuumDecay/reference/v1.0/*`: frozen legacy regression references.
  * All are kept, moved to `legacy/` with a data dictionary. None is deleted.
* Delete the AI-prompt documents (`docs/prompt_for_claude.md`, `docs/critique_prompt_for_claude.md`,
  `docs/reformulate_prompt_for_claude.md`) and the AI-agent configuration `docs/AGENTS.md`.
* Delete `docs/README` (the Bertini examples README), after confirming that nothing in the
  repository uses Bertini.
* Move the other predecessor planning/review documents to `docs/archive/` with an index
  marking them historical and superseded.

## 10. Licence (owner decision 6)

Untouched. The repository has no licence file at the top level. `SMVacuumDecay/LICENSE` is
empty. `LICENSES/` holds Bertini and GMP/MPFR licence texts for software the repository does not
use. These files are licensing metadata, so they are **left in place unchanged**. The ambiguity
is documented for the owner to decide.
