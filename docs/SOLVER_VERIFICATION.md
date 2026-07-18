# Solver Verification: Code-Level Audit

**Role:** Hostile PRD referee  
**Target:** `solver_numerical.cpp`  
**Supporting:** `paper/main.tex`, `solver_analytical.cpp`  
**Date:** 2026-06-16

---

## 1. Is the bounce profile obtained by solving the Euclidean bounce ODE with boundary conditions?

**Answer: No.**

The code never calls an ODE solver for the bounce equation. The bounce profile is
computed from a **closed-form ansatz** at `solver_numerical.cpp:327`:

```cpp
double phi_x  = prefactor * 2.0 / (R * (e2x + 1.0));
```

where (`solver_numerical.cpp:353–354`):

```cpp
double prefactor = std::sqrt(2.0 / abs_lam);
```

This is the Fubini-Lipatov instanton profile:

$$\phi(r) = \sqrt{\frac{2}{|\lambda_R|}} \; \frac{2R}{r^2 + R^2}$$

The transformation `u = (e^{2x} - 1) / (e^{2x} + 1)` maps to the compactified
integration variable; the profile at `solver_numerical.cpp:327` is algebraically
identical.

**There is no ODE solve.** The Euclidean bounce equation (`main.tex` Eq. 19,
`solver_numerical.cpp` lines 219–223 in pseudocode references):

$$ \frac{d^2\phi_b}{dr^2} + \frac{3}{r}\frac{d\phi_b}{dr} = V'_{\rm eff}(\phi_b) $$

is **never discretized, never integrated, never solved**. The code contains no
shooting method, no relaxation, no collocation — nothing that would produce a
numerical solution to this equation.

---

## 2. Is a Fubini-Lipatov ansatz assumed with only the radius optimized?

**Answer: Yes — exactly and exclusively.**

The evidence:

**(a) Kinetic term is closed-form, not integrated** (`solver_numerical.cpp:355`):

```cpp
double kinetic_term = (16.0 * pi * pi) / (3.0 * abs_lam);
```

This is the **analytical** kinetic energy of the Fubini-Lipatov profile for a
pure quartic potential:

$$S_{\rm kin} = 2\pi^2 \int_0^\infty dr\, r^3 \frac{1}{2} \left(\frac{d\phi}{dr}\right)^2
= \frac{16\pi^2}{3|\lambda_R|}.$$

The kinetic contribution is **never evaluated numerically**. It is inserted by
hand from the known closed-form result for the Fubini-Lipatov profile with
**constant** $\lambda = \lambda_R$.

**(b) Only the potential integral is numerical** (`solver_numerical.cpp:359–368`):

```cpp
double potential_integral = 0.0;
for (int i = 0; i <= N; ++i) {
    double u = -1.0 + i * du;
    ...
    double f = integrand_u(rge, u, R, mu_inst, prefactor);
    ...
    potential_integral += weight * f * du;
}
return kinetic_term + potential_integral;
```

This integrates $V_{\rm eff}(\phi(r))$ **using the Fubini-Lipatov profile**
$\phi(r)$ (line 327). The potential integral is the only numerically evaluated
component of the bounce action.

**(c) The only optimized degree of freedom is $R$** (`solver_numerical.cpp:372–415`):

The Golden-Section Search (GSS) varies a **single parameter**: the bubble radius $R$.
There is no shape optimization — no variation of the profile function, no
waveform relaxation, no basis-function expansion.

**(d) The function `rk4_adaptive_step` is defined (lines 282–304) but never called.**
The main loop (`solver_numerical.cpp:446–470`) calls only `rk4_single_step`
(fixed dt=0.1) for the RGE integration, not for the bounce.

---

## 3. Which equations in the manuscript would become inaccurate?

| Manuscript equation | Problem |
|---|---|
| **Eq. (19)** — The bounce EOM $\phi_b'' + (3/r)\phi_b' = V_{\rm eff}'(\phi_b)$ | **Never solved.** The paper presents this ODE as the defining equation for $\phi_b(r)$, then never uses it. |
| **Eq. (20)** — $S_E = 2\pi^2 \int r^3 [\frac{1}{2}\phi_b'^2 + V_{\rm eff}(\phi_b)]\, dr$ | **Partially evaluated.** The kinetic term uses a closed-form expression valid only for constant $\lambda$. The integrand uses a Fubini-Lipatov profile, not a solution of Eq. (19). |
| **Eqs. (22)–(24)** — Numerical error budget | **Groundless.** The $10^{-14}$ RK4 truncation error (line 489) applies only to the RGE solver, not the bounce ODE (which is never integrated). The $10^{-11}$ quadrature error (line 490) is the error of Simpson's rule on the *Fubini-Lipatov* integrand, not a true bounce. |
| **§IV.D Eq.-of-motion residuals** | **Circular.** The "residual" checks how well the Fubini-Lipatov profile satisfies the bounce ODE, but the profile was never designed to solve it. A non-zero residual does **not** indicate numerical error — it indicates a wrong profile. |

---

## 4. What claims in the paper must be softened or rewritten?

### Critical (must change)

1. **Abstract, line 36–43:**  
   "By deploying an **exact numerical framework** integrating three-loop
   renormalization group equations with a 4th-order Runge-Kutta solver...
   we achieve numerical errors below $10^{-11}$ in relative action"

   **Problem:** The framework is not "exact numerical." It is a
   Fubini-Lipatov ansatz with a numerically evaluated potential correction.
   The $10^{-11}$ error figure is misleading because it excludes the
   **systematic error of the profile ansatz itself**, which Appendix D admits
   can be $\mathcal{O}(1)$ or divergent (negative action).

2. **Abstract, line 52:**  
   "Our results demonstrate that **exact numerical bounce solvers** are strictly
   necessary for precision vacuum stability bounds within the $1\sigma$
   experimental uncertainty of the top quark mass."

   **Problem:** This paper does NOT use an exact numerical bounce solver.
   It uses a semi-analytical ansatz. The statement that such solvers "are
   strictly necessary" is a logical conclusion that may be correct, but the
   paper does **not** demonstrate it, because it does not deploy one.

3. **§III (Numerical Methodology), entire section:**  
   The section title and content imply a full numerical solution. The words
   "bubble profile" (line 271), "numerical evaluation of the Euclidean action
   integral for a given bubble profile" (line 271–272), and "optimization of
   the bubble radius $R$" (line 272) collectively create the impression of a
   general-purpose solver. The fact that the profile is a fixed ansatz is
   disclosed only in Appendix D ($\S$IV.D, line 749).

### High (should change)

4. **§III.B, line 296–312 (Quadrature section):**  
   Describes Simpson's rule for the action integral but never states that the
   integrand uses the Fubini-Lipatov profile. The substitution $u = \tanh(r/R)$
   is described (line 301) but its purpose (to compactify the closed-form
   profile) is not connected to the ansatz.

5. **§IV.D, Eq.-of-motion residuals, lines 444–468:**  
   The claim that "numerical bounce profiles tightly satisfy their defining
   ODE throughout the parameter space" (lines 452–454) is misleading. The
   Fubini-Lipatov profile is an **approximate** solution of the bounce ODE —
   exact only for a pure quartic with constant $\lambda$. The residual
   distribution is a measure of **how wrong the ansatz is**, not a validation
   of the numerics.

### Should change (moderate)

6. **Table I (error budget), lines 486–498:**  
   Lists RK4 truncation, quadrature, and GSS errors without any entry for
   **ansatz systematic error**. The dominant error source — the Fubini-Lipatov
   approximation — is completely absent from the budget table.

7. **§III.A line 284–294 (RK4 description):**  
   "The adaptive step-size controller enforces that the local truncation error
   at each step satisfies $|\delta\phi| < \epsilon_{\rm RK4}$" (lines 290–293).

   **Problem:** `rk4_adaptive_step` is defined but **never called**.
   The RGE integration uses fixed dt=0.1 (line 433). The paper's description
   of "adaptive step sizing" is factually incorrect about the code.

---

## 5. Is it still valid to describe the method as "exact numerical bounce action computation"?

**No.**

The term "exact" in computational physics has a specific meaning: the numerical
solution converges to the true solution of the governing equations in the limit
of vanishing discretization error. This paper's method does **not** converge to
the true bounce solution in any limit, because:

1. The profile is **constrained** to the Fubini-Lipatov family
   $\phi(r) = A \cdot 2R/(r^2 + R^2)$ with $A = \sqrt{2/|\lambda_R|}$ fixed.
2. Refining the Simpson grid ($N \to \infty$) or the GSS tolerance
   ($\tau \to 0$) does **not** change the profile shape — it only evaluates
   the action integral more precisely for a predetermined, possibly wrong,
   profile.
3. The paper's own Appendix D (line 749) admits the profile is not a saddle
   point when $\lambda$ runs, which is precisely the regime of interest.

The method is best described as a **semi-analytical bounce action estimator**
using a **Fubini-Lipatov profile ansatz** with **numerical RGE running** and
**numerical potential integration**.

---

## 6. What terminology would be technically correct?

| Current | Correction |
|---------|-----------|
| "exact numerical bounce solver" | "semi-analytical bounce action computation based on the Fubini-Lipatov ansatz" |
| "exact numerical framework" | "hybrid numerical-analytical framework" |
| "exact bounce action computation" | "bounce action computation using the conformal profile with running-coupling corrections" |
| "numerically exact bounce profiles" | "Fubini-Lipatov profile with Simpson-evaluated potential integral" |
| "bounce profile" (unqualified) | "Fubini-Lipatov profile (see Appendix D)" on first use |
| "numerical error budget" | "numerical error budget for the RGE integration, quadrature, and GSS (excluding ansatz systematics)" |
| "adaptive step-sizing" (for RK4) | "fixed-step RK4 with discussed-in-text tolerance" (or implement adaptive stepping) |

---

## Appendix: Key Code Quotations

| What | File:Line | Code |
|------|-----------|------|
| Fubini-Lipatov profile | `solver_numerical.cpp:327` | `double phi_x = prefactor * 2.0 / (R * (e2x + 1.0));` |
| Prefactor = $\sqrt{2/|\lambda_R|}$ | `solver_numerical.cpp:353` | `double prefactor = std::sqrt(2.0 / abs_lam);` |
| Kinetic term = closed form | `solver_numerical.cpp:355` | `double kinetic_term = (16.0 * pi * pi) / (3.0 * abs_lam);` |
| Only potential integral is numerical | `solver_numerical.cpp:359–368` | Simpson loop; only $V_{\rm eff}$ terms appear |
| Fixed RGE step (no adaptive) | `solver_numerical.cpp:433` | `double dt = 0.1;` |
| RGE integration uses fixed step | `solver_numerical.cpp:449` | `y = rk4_single_step(y, t, dt);` (NOT `rk4_adaptive_step`) |
| Adaptive function defined but unused | `solver_numerical.cpp:282–304` | `bool rk4_adaptive_step(Params& y, double& t, double& dt)` |
| Only $R$ is optimized | `solver_numerical.cpp:372–415` | GSS varies `R` only; no shape degrees of freedom |
| **No bounce ODE solve anywhere** | entire file | No shooting, relaxation, or collocation of Eq. (19) |
