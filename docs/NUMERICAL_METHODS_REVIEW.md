# Numerical Methods Review

**Project:** SMVac — Naman Goyal (namangoyal31), IISc Bangalore
**Status:** Stage B review of every numerical algorithm in the original implementation
(commit `a7c070a`). Each claim is backed by a script in `audit/` that can be re-run; see
`audit/README.md`. Physics issues are discussed in `docs/SCIENTIFIC_METHODOLOGY.md`. This
document separates *numerical* error from *method/physics* error.

> **Paths** refer to the reviewed revision `a7c070a`. The original code, data and figures are now in `legacy/`,
> and its working notes from `docs/` are in `docs/archive/`. The replacement methods chosen afterwards are
> documented in `docs/NUMERICAL_METHOD_COMPARISON.md` and validated in `docs/VALIDATION_REPORT.md`.

---

## 1. Inventory of algorithms

| Algorithm | Location | Purpose |
|---|---|---|
| Classical RK4, fixed step Δt = 0.1 in $t=\ln\mu^2$ (≈780 steps from $M_t$ to $M_{\rm Pl}$) | `rk4_single_step`, loop in `classify_stability` | RG evolution |
| Table lookup: linear interpolation of $\lambda$ only, other couplings piecewise constant | `RGEHelper::get_params` | couplings at arbitrary $\mu$ |
| Minimum and zero crossing of $\lambda_{\rm eff}$ taken on grid nodes | `classify_stability` | $\lambda_{\rm eff,min}$, instability scale |
| Composite Simpson, $N=2048$, on $u\in(-1,1)$ with $r=R[(1+u)/(1-u)]^{2}$ | `evaluate_action_at_R` | potential integral of the FL profile |
| Golden-section search over $\ln R\in[\ln R_*-8,\ln R_*+8]$, tol $10^{-13}$ | `find_minimum_action` | minimise `S_exact`(R) |
| Unused adaptive RK4 (step doubling, tol $10^{-13}$) | `SMVacuumDecay/numerics/RK4.hpp` | not called anywhere |

## 2. RG integration — adequate as is [verified]

The RGEs are smooth and only mildly stiff: the coefficients are $\mathcal O(10^{-2})$ per unit $t$.
At the central point $(125.15,173.34)$ (`audit/exp_defect_impact.py`):

| Δt | $\lambda(M_{\rm Pl})$, RK4 − DOP853 (rtol $10^{-12}$) |
|---|---|
| 0.2 | $+1.37\times10^{-12}$ |
| 0.1 (legacy) | $+8.5\times10^{-14}$ |
| 0.05 | $+4.9\times10^{-15}$ |

The error ratio of 16 per halving confirms fourth order. The legacy integrator is accurate to
$\sim10^{-13}$, about $10^{9}$ below the physics uncertainty ($\pm3\times10^{-4}$ on $\lambda(M_t)$).
Adaptive or stiff solvers are **not** needed for accuracy. One small defect remains: the
perturbativity guard tests only $\lambda$ and $y_t$.

*Recommendation:* keep fixed-step RK4 (simple, reproducible, demonstrably converged). Store
$\beta(t)$ at every node, so the couplings can be interpolated with cubic Hermite polynomials
(next section).

## 3. Coupling interpolation — introduces roughness [verified]

`get_params` interpolates $\lambda$ linearly but returns the other six couplings from the left
node. The $\lambda_{\rm eff}(h)$ fed to the bounce integrand is therefore discontinuous at every
node, with jumps of $\mathcal O(|dy_t/dt|\,\Delta t\,\partial\lambda_{\rm eff}/\partial y_t)$, up to
$\sim10^{-4}$ at low scales. Measured effect on the `S_exact` objective
(`audit/exp_action_landscape_zoom.py`): within $|\Delta\ln R|\le0.02$ of the optimum, the first
difference of $S(R)$ changes sign **509–636 times out of 2000** samples. The roughness is
$1.5\times10^{-6}$ to $5\times10^{-5}$ relative.

*Recommendation:* cubic Hermite interpolation of **all** couplings, using the stored exact
derivatives. This is $C^1$ and $\mathcal O(\Delta t^4)$ accurate, with no extra RGE evaluations.
`audit/independent.py` uses a DOP853 dense output for the same purpose.

## 4. Grid-node extrema — adequate, but replaceable

The grid minimum of $\lambda_{\rm eff}$ agrees with a continuous minimisation to about
$3\times10^{-11}$ (Δt = 0.05–0.2 comparison), because the minimum is flat. The instability scale
has a resolution of Δln μ = 0.05, i.e. 5% in μ.

*Recommendation:* use Brent root-finding and minimisation on the smooth interpolant. The cost is
negligible, and it removes grid dependence from reported scales.

## 5. Quadrature of the FL potential integral

The algebraic map $r=R[(1+u)/(1-u)]^{\alpha/2}$ with α = 4 makes the integrand vanish like
$(1\mp u)^7$ at both ends, so the endpoints are harmless. With a smooth integrand, Simpson at
$N=2048$ would indeed be very accurate. The integrand, however, inherits the node
discontinuities of §3 (and the frozen-coupling kinks at $v$ and $M_{\rm Pl}$), so the formal
$\mathcal O(h^4)$ convergence claimed in `docs/error_analysis_report.md` does not apply. More
importantly, this quadrature evaluates an unsound objective (SCIENTIFIC_METHODOLOGY §6.3), so
improving it is pointless.

*Recommendation:* retire it together with `S_exact`. The replacement bounds need
$\bar\lambda(h)=6\int_0^1w(1-w)\lambda(hw)\,dw$. In $x=\ln w$ this is a smooth one-dimensional
integral, done with adaptive Gauss–Kronrod (QUADPACK `qags` in the audit; $10^{-13}$ requested).

## 6. Golden-section search — mis-specified [verified]

* **Tolerance.** $10^{-13}$ in $\ln R$ is meaningless. With noise of relative size ε on a
  quadratic minimum, the attainable resolution in the argument is $\sim\sqrt\varepsilon$, and
  here the objective is rough at $10^{-6}$–$10^{-5}$ (§3). The 70-odd iterations beyond that
  point only random-walk on the noise.
* **Platform dependence.** The same source compiled for x87 (80-bit) and for SSE2 (64-bit)
  arithmetic gives `S_exact` differing at $\sim6\times10^{-7}$ relative at
  $(M_h,M_t)=(125,173)$ (2433.39597 vs 2433.39753) and $2\times10^{-7}$ at (115, 180). The SSE2
  build agrees exactly with the Python port. The tests assert agreement to $10^{-12}$, which
  holds only for the benchmark points and compiler that generated the references.
* **Non-unimodality and discontinuities.** The objective has a $10^{100}$ plateau (λ_c ≥ 0) and a
  jump at $T=5\times10^5$, where the code returns $T$ instead of $T+V$. It can also be unbounded
  below (§6.3 of the methodology). The search converged onto the jump at (5, 96) and onto a cliff
  at (5, 105) (`audit/results/action_landscape.png`). At (5, 96), 696 of 801 bracket points are the
  sentinel.
* **Success criteria.** The routine always returns a number. Failure (sentinel, early exit,
  boundary minimum) is never reported, and the result is fed straight into the classification.

*Recommendation:* retire this search with its objective. For $\min_h\bar\lambda(h)$, which is
smooth and unimodal near the minimum, use a coarse scan followed by Brent, with tolerance
$\sim10^{-8}$ in $\ln h$, and report failures explicitly.

## 7. Bounce computation — proposed methods

| Method | What it gives | Cost | Validation |
|---|---|---|---|
| FL estimate $8\pi^2/(3\lvert\lambda_{\min}\rvert)$ | rigorous lower bound (Pohozaev + sharp Sobolev); leading-order estimate | one minimisation | constant-λ limit |
| CGM variational estimate $8\pi^2/(3\lvert\min_h\bar\lambda\rvert)$ | upper estimate (rigorous only if the CGM infimum is attained) | one 1-D minimisation of a quadrature | equals lower bound for constant λ; checked to lie above the numerical bounce |
| Shooting in $(\tau,f)$ with bisection on $\ln\phi(0)$ | numerical solution of the tree-level O(4) bounce equation | ≈10–100 ODE solves | Derrick and dilatation identities; tolerance refinement; lies between the bounds |

For the pure-SM potential the gap between the lower bound and the upper estimate is below 0.1% at all benchmark points
(table in SCIENTIFIC_METHODOLOGY §6.6). That is far below the theory systematics of a tree-level
treatment. The bounds therefore suffice for the phase diagram, and shooting serves as an
independent check.

*Why not a relaxation or BVP solver.* The bounce of a nearly scale-invariant potential has a
quasi-zero mode (dilatation), which makes collocation and relaxation Jacobians nearly singular.
Shooting in log variables handles this naturally, because the mode becomes the selection of
$\phi(0)$ by the dilatation identity.

## 8. Claimed versus actual numerical accuracy

| Claim (source) | Status |
|---|---|
| "adaptive RK4, ε = 1e-13" (`docs/error_analysis_report.md`, `SMVacuumDecay/README.md`) | false: fixed step is used. The fixed step is nevertheless accurate to $10^{-13}$. |
| "quadrature error < 1e-11" | not applicable: the integrand is discontinuous, and the objective is unsound |
| "total numerical error ≲ 0.1%, classification exact > 0.5 GeV from boundaries" | unsupported: method error is $\mathcal O(1)$ and the sign can be wrong; physics defects shift boundaries by ≈0.3 GeV, and the threshold is off by 45–65 units |
| "OpenMP parallelisation" (README) | false: no OpenMP pragmas. With MinGW.org GCC 6.3, `-fopenmp` fails to link (missing `-lpthread`). |
| "Exact evaluation of the Euclidean action" (`SMVacuumDecay/README.md`) | false |

## 9. Prioritised improvement plan (impact × risk)

| Priority | Change | Scientific impact | Risk | Test |
|---|---|---|---|---|
| P1 | Retire `S_exact`; classify with the FL bounds (+ shooting validation) | removes O(1) errors and sign failures | low | constant-λ limit, bounds bracket shooting, identities |
| P2 | Correct $g_2$ (and $g_1$, $y_t$ slopes) to Buttazzo et al. eqs. (55)–(60); make inputs configurable | Δλ(M_Pl) = −1.5e-3 | low | reproduce eq. (61) |
| P3 | Fix the three RGE transcription errors | Δλ(M_Pl) = −3.5e-4 | low | power-counting test; eq. (61) |
| P4 | Fix the Z coupling in λ_eff | ΔS ≈ 3.5% | low | analytic CW check |
| P5 | Literature metastability criterion (bounce scale, $0.15/H_0^4$) | boundary shift | low | analytic formula test |
| P6 | Remove $\phi^6$ from the SM baseline (optional explicit BSM parameter) | removes the 2.4% bias | decision | toggled comparison |
| P7 | Restrict the primary scan to the validity domain of the matching | removes unsupported regions | decision | — |
| P8 | Hermite interpolation; Brent extrema; explicit failure statuses | smoothness, honest errors | low | convergence tests |
| P9 | Build system (CMake), unit/physics/regression tests, config files, metadata in outputs | reproducibility | low | CTest |
