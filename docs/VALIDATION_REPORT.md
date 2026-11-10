# Validation Report (Stage D)

**Project:** SMVac — Naman Goyal (namangoyal31), IISc Bangalore

This report demonstrates the reliability of the numerical reference and of the inputs it depends on. Every
number can be regenerated with the scripts named in each section; quoted scan statistics come from
`results/summary.json` (`validation/analyze_scans.py`).

Throughout, *numerical reference* means a numerical solution with an estimated uncertainty. Nothing here is
called exact. *Rigorous* is used only for the lower bound $S\ge S_{\rm FL}$.

---

## 1. Physics inputs

| Check | Reference | Result | Where |
|---|---|---|---|
| Matching at $M_t$ (central values) | Buttazzo et al. (2013) eqs. 55–60 | identical by construction (unit-tested to $10^{-12}$) | `python/tests`, `src/matching.cpp` |
| Loop-order power counting of all β-functions | degree $n+1$ at $n$ loops (literature independent) | all pieces consistent to $10^{-13}$; the same test flags exactly the two 3-loop errors of the legacy code | `tests/test_beta_functions.cpp`, `audit/check_power_counting.py` |
| One-loop coefficients, four-loop QCD sign | standard SM expressions; $-\beta_3(n_f=6)=-2472.28$ (van Ritbergen et al. 1997) | pass | `tests/test_beta_functions.cpp` |
| Couplings at $M_{\rm Pl}$ | Buttazzo et al. eq. (61): $g_1=0.6154$, $g_2=0.5055$, $g_3=0.4873$, $y_t=0.3825$, $\lambda=-0.0143$ | C++: 0.615367, 0.505459, 0.487315, 0.382499, −0.0143688 (within the quoted precision) | `tests/test_rge.cpp` |
| Independent transcription | Python reference (legacy expressions + corrections, DOP853) vs C++ (fresh transcription, RK4) | $\lambda(M_{\rm Pl})$ agrees to $<10^{-10}$; $S_{\rm FL}$ to $10^{-13}$; $S_{\rm CGM}$ to $10^{-13}$ | `tests/test_rge.cpp`, `tests/test_analytic.cpp` |
| Coleman–Weinberg $\lambda_{\rm eff}$ | explicit $\sum_i n_iM_i^4(\ln M_i^2/\mu^2-c_i)/64\pi^2$ with $M_Z^2=(g_2^2+g_Y^2)h^2/4$ | agrees to $10^{-12}$; the legacy $g_1^2$ variant is rejected | `tests/test_potential.cpp` |

## 2. RG integration and interpolation

* Classical RK4 converges at fourth order: the error ratio is ≈16 per step halving (`tests/test_rge.cpp`;
  NUMERICAL_METHODS_REVIEW §2). With step 0.025 the integration error in $\lambda(M_{\rm Pl})$ is $<10^{-15}$.
* Cubic-Hermite dense output is accurate to $<10^{-9}$ at mid-nodes even for step 0.1, and to $\sim10^{-13}$ at
  the production step.
* The Chebyshev representation of $\lambda_P(x)$ is truncated where the tail falls below $10^{-14}$. The
  sensitivity of $S$ to the RG tolerance and to the Chebyshev degree is in §4.

## 3. Bounce solver

### 3.1 Analytic limit

Constant $\lambda$: the solver reproduces $8\pi^2/(3|\lambda|)$ to $10^{-10}$ (`tests/test_bounce.cpp`). This
problem is exactly degenerate (every translate is a solution); the residual-based Newton stop handles it.

### 3.2 Cross-method agreement

On every benchmark where they converge, the C++ M4 solver and the independent Python solvers M1 (shooting,
physical variables), M2 (shooting, scale-covariant variables) and M4 (Chebyshev) agree to
$\lesssim10^{-11}$ relative (docs/NUMERICAL_METHOD_COMPARISON.md). The C++ tests pin:

* parabola $q=250$: $2107.0357973921$, agreement $<10^{-10}$;
* parabola $q=2.5$: $225400.864807$, agreement $4\times10^{-13}$;
* SM central point (`cw1`): $1856.01513783$, agreement $2\times10^{-12}$.

### 3.3 Identities and self-convergence

For every accepted solution:

* Derrick and dilatation residuals are $<10^{-9}$ (standard tier) or $<10^{-5}$ (reduced-accuracy tier);
* the two-resolution uncertainty $\delta S$ is reported per point.

Scan statistics are in §5.

### 3.4 Bounds

The rigorous lower bound must hold for every solution, and the variational estimate should lie above it:
$S_{\rm FL}\le S_{\rm num}\le S_{\rm CGM}$. Violations over the scans are counted in §5.

## 4. Convergence, precision and residual checks

The checks below do not rely on agreement between solvers.

### 4.1 Self-convergence in the discretisation

Source: `validation/results/convergence.csv` (`validation/convergence_study.py`). The table gives the relative
deviation from the shooting consensus at rtol $=10^{-13}$, and from M4 at $N\ge192$ where that is more accurate.

| case | M4, $N=256$ | M4, $N=384$ | M4, $N=768$ | M4 plateau | M2 / M1, rtol $10^{-10}$ | M2 / M1, rtol $10^{-13}$ |
|---|---|---|---|---|---|---|
| parabola $q=250$ (thick wall) | 5e-12 | 5e-12 | 5e-12 | constant for $N\ge192$ | 6e-14 / 5e-12 | 6e-12 / 6e-12 |
| parabola $q=2.5$ (thin-wall onset) | 2e-12 | 2e-12 | 2e-12 | constant for $N\ge192$ | 1e-12 / 8e-12 | 6e-12 / 6e-12 |
| parabola $q=0.1$ (thin wall) | 6e-05 | 9e-08 | 8e-13 | reached at $N=768$ | 2e-12 / 9e-12 | 1e-12 / 1e-12 |
| SM `cw1` (125.15, 173.34) | 5e-12 | 5e-12 | 5e-12 | constant for $N\ge192$ | 2e-11 / 4e-12 | 2e-12 / 2e-12 |
| SM `cw1` (125.15, 171.149), near-critical | 3e-12 | 8e-12 | 6e-12 | constant for $N\ge192$ | 7e-12 / 6e-12 | 3e-12 / 3e-12 |

* **The plateau is not M4's error.** The M4 deviation is constant in $N$ once converged, so it measures the
  accuracy of the shooting reference itself (≈ $5\times10^{-12}$, limited by quadrature of its dense output). M4's
  own self-convergence is at the round-off level (§4.2).
* **Thin walls need more nodes.** At $q=0.1$, fixed $N=256$ is insufficient ($6\times10^{-5}$). The production
  solver detects this through the identity residuals and doubles $N$; this is the *adaptive* rule of
  `src/bounce.cpp`.
* **The shooting methods** converge roughly linearly in the integrator tolerance down to that plateau.

### 4.2 Floating-point precision

Source: `validation/precision_check.cpp`, output `validation/results/precision_check.csv`. The same collocation
discretisation is solved in `double` and in `long double` (64-bit mantissa), with the potential evaluated
natively in the wider type.

| case | $N$ | $(S_{\rm long\,double}-S_{\rm double})/S$ |
|---|---|---|
| parabola $q=250$ | 192 / 256 / 384 | 1.3e-15 / 1.9e-15 / 0 |
| parabola $q=2.5$ | 192 / 256 / 384 | 3.6e-15 / 1.9e-15 / 2.6e-16 |
| SM `cw1` central | 192 / 256 / 384 | 2.6e-15 / 1.2e-16 / 1.8e-15 |

The round-off contribution of the double-precision solver is $\lesssim4\times10^{-15}$ relative. The Derrick
residual at the SM point (~$6\times10^{-13}$) is the same in both precisions, so it is not a round-off effect.
The values are from the CMake target `precision_check` (SSE2 double, `-O3`). A build without SSE2 on 32-bit x86
evaluates the double path in x87 extended precision; it gives differences of the same order ($\le1.3\times10^{-15}$).

### 4.3 Domain truncation, off-node residual, quadrature rule, input round-off

Source: `validation/results/discretization_checks.txt` (`validation/discretization_checks.py`).

| check | parabola $q=2.5$ | SM central | SM near-critical (S = 1.41e5) |
|---|---|---|---|
| half-width $L=10\to25$ at fixed node density: change of $S$ | $<10^{-15}$ | $<10^{-14}$ | $<10^{-14}$ |
| off-node ODE residual, max / rms (relative to max f) | 1.1e-8 / 9e-11 | 2.8e-8 / 2.7e-10 | 2.1e-8 / 1.6e-10 |
| Simpson (20001 points) vs Clenshaw–Curtis action | 4e-14 | 4e-15 | 1.5e-12 |
| 1-ulp random perturbation of the $\lambda(x)$ coefficients (spread of $S$) | – | 3.6e-15 | 2.7e-13 |

**Boundary conditions.** The Robin conditions $f'=\pm f$ hold exactly only asymptotically. The neglected term
$f^3G$ is suppressed by $e^{-2|\tau-\tau_c|}$ relative to the linear one. The domain check (first row) shows that
the truncation at half-width 15 changes $S$ by less than $10^{-14}$. The two integral identities (Derrick,
dilatation) are monitored at every point (§3.3). They are not imposed by the discretisation.

The collocation residual vanishes at the nodes by construction. Evaluated on a 20001-point grid from the
Chebyshev interpolant, its maximum sits at the round-off floor expected for second derivatives of a 384-term series,
$\varepsilon N^4/L^2\approx2\times10^{-8}$. The rms is two orders of magnitude smaller.

### 4.4 Physics-input representation

For both SM test points, changing the RG tolerance ($10^{-13}\to10^{-10}$) changes $S$ by $\le5\times10^{-11}$.
Changing the Chebyshev degree of $\lambda(x)$ (128, 256, 512) changes it by $\le10^{-12}$
(`validation/results/convergence.csv`).

### 4.5 Uncertainty assigned to the reference

Combining §4.1–4.4, the numerical uncertainty of $S_{\rm num}$ is:

* about $10^{-14}$ relative at generic points;
* about $10^{-12}$ near criticality, where round-off in the inputs is amplified;
* the per-point estimate $\delta S$ for reduced-accuracy points, which is typically $10^{-6}$ relative.

The production $\delta S$ (difference between two resolutions plus identity residuals times $S$) is a
conservative upper estimate. At all standard-tier points of the controlled scans it is $\le6.7\times10^{-10}$
relative (`cw1`) and $\le1.1\times10^{-9}$ (`tree`).
These uncertainties are far below every physics effect discussed in this report.

## 5. Scan-level checks (controlled region)

Sources:

* `results/controlled_cw1.csv`, `results/controlled_tree.csv`, `results/controlled_cw1_fallback.csv` and
  `results/summary.json`, from `python/smvac/run_scan.py`, `validation/fallback_shooting.py` and
  `validation/analyze_scans.py`. For their provenance (`git_commit = bb4b5b7…-dirty`), see README §8. In short:
  * the scan binary was built from `bb4b5b7` plus one uncommitted change, `max_nodes = 1024`, later committed in
    `ad63002`;
  * other source edits made during the scans did not reach that binary;
  * a 233-point recomputation from `ad63002` is supporting evidence, not a full regeneration;
* `validation/results/near_critical_*.csv` and `near_critical_summary.json`, from `validation/near_critical.py`
  and `validation/analyze_near_critical.py`.

Grid: $M_h\in[120,130]$, $M_t\in[166,182]$ GeV in steps of 0.1 GeV (16,261 points per prescription).

### 5.1 Solver performance and per-point checks

| | `cw1` | `tree` |
|---|---|---|
| stable / metastable / unstable (collocation) | 5117 / 7052 / 3555 | 5067 / 7519 / 3675 |
| no finite bounce, $\lambda_{\min}$ on the frozen boundary at $M_{\rm Pl}$ | 177 | 0 |
| collocation failures | 360 (all resolved by the fallback, §5.2) | 0 |
| reduced-accuracy tier (identities $<10^{-5}$; Planck kink) | 164 | 0 |
| flagged Planck-sensitive | 351 | 0 |
| max $\delta S/S$, standard tier, not Planck-sensitive | $6.7\times10^{-10}$ | $1.1\times10^{-9}$ |
| max identity residual, same points | $5.0\times10^{-10}$ | $8.3\times10^{-10}$ |
| violations of $S_{\rm FL}\le S_{\rm num}\le S_{\rm CGM}$ | 0 | 0 |

* **Planck sensitivity is a `cw1` effect.** In `tree` the minimum of $\lambda$ lies below $10^{18.64}$ GeV
  throughout the region, so no `tree` point needs the fallback or the reduced tier. In `cw1` it reaches
  $M_{\rm Pl}$ for $M_h\gtrsim127$ GeV.
* **Prescription systematic.** `tree` and `cw1` differ in $S$ by 2.6–5.0% at 90% of the common bounce points
  (median 2.7%). The difference grows near criticality, because $M_t^{\rm crit}$ itself differs by 0.05 GeV. The
  two prescriptions classify 707 of 16,261 points differently.

### 5.2 Collocation failures: shooting fallback and the Planck plateau

All 360 `cw1` collocation failures lie at $M_h\ge127$ GeV, where the minimum of $\lambda_{\rm eff}$ approaches
$M_{\rm Pl}$ and the bounce meets the frozen-coupling kink. Each was re-solved by the two independent shooting
solvers M2 and M1. Acceptance required convergence of both, $|S_{\rm M2}-S_{\rm M1}|<10^{-8}S$ and
$|{\rm Derrick}|<10^{-8}$.

* **78 points are metastable.** M2 and M1 agree to $\le5.3\times10^{-10}$ and $|{\rm Derrick}|\le5.0\times10^{-10}$.
  All 78 are Planck-sensitive, with $\phi(0)=10^{19.15}$–$10^{20.74}$ GeV.
* **282 points have no finite bounce.** The minimum of $\lambda_{\rm eff}$ lies at $10^{18.87}$–$10^{19.08}$ GeV,
  just below $M_{\rm Pl}=10^{19.09}$ GeV, and $\lambda$ rises by at most $7.1\times10^{-6}$ before the plateau.
  Both shooting solvers find only undershoots over the whole admissible interval, and the CGM optimum of the FL
  family escapes to arbitrarily large field. The infimum of the action is approached only by bubbles that escape
  into the plateau. It equals $8\pi^2/(3|\lambda(M_{\rm Pl})|)\le1.024\,S_{\rm FL}$. This is the same mechanism
  as for the 177 points whose minimum sits exactly on the boundary (METHODOLOGY §6.6, last row).

**This is a property of the frozen-coupling model, not a numerical failure.** The two are told apart by explicit
criteria:

* a *no-finite-bounce* point needs every shot of both solvers to undershoot, and the variational optimum to
  run to the edge of field space;
* a *solver failure* is any other non-convergence or disagreement.

After the fallback, the controlled `cw1` scan contains **no unresolved points**:

| status | stable | metastable | unstable | no finite bounce (Planck plateau) |
|---|---|---|---|---|
| points | 5117 | 7130 | 3555 | 459 = 177 + 282 |

**Classification of the plateau points.** All 459 have $S_{\rm FL}\ge1301$. The lower bound is a statement about
*critical points*, not about arbitrary field configurations, on which the action is unbounded below (METHODOLOGY
§6.2). It enters here as follows.

* *Mathematical step (rigorous for the frozen-coupling potential).* For any configuration with $V_{\rm int}<0$,
  the scale-maximised action satisfies $\max_aS(\phi_a)=T^2/(4|V_{\rm int}|)\ge S_{\rm FL}$, by the sharp Sobolev
  inequality.
* Every path from the false vacuum to a configuration with negative action crosses the Pohozaev manifold
  $T+2V_{\rm int}=0$, where $S=T^2/(4|V_{\rm int}|)$. The mountain-pass barrier therefore satisfies
  $S_{\rm barrier}\ge S_{\rm FL}$. For these points, the infimum approached by escaping bubbles,
  $8\pi^2/(3|\lambda(M_{\rm Pl})|)$, is itself $\ge S_{\rm FL}$.
* *Physical assumptions (not rigorous).* The semiclassical formula
  $\wp=0.15\,\Lambda_B^4e^{-S}/H_0^4$ is used with the frozen-coupling potential and with $\Lambda_B\le M_{\rm Pl}$.
  Then $\wp<1$ whenever $S>S_{\rm crit}=\ln0.15+4\ln(M_{\rm Pl}/H_0)=559.3$.

Under these assumptions the 459 points are metastable, because $S_{\rm FL}\ge1301>S_{\rm crit}$. This conclusion
depends on freezing the couplings above $M_{\rm Pl}$ and on the semiclassical, tree-level, non-gravitational
treatment. Their lifetime is set entirely by the treatment above $M_{\rm Pl}$, so no SM prediction is made for
them.

### 5.3 Near-critical behaviour

**Model family.** $\lambda=\lambda''[\tfrac12(x-x_0)^2-q]$ with 17 values of $q$ from 1000 down to 0.05. The
reference is the Python M4 with adaptive $N$ (256 for $q\ge1$, up to 1024 for $q\le0.1$), cross-checked by M2 to
$\le1.1\times10^{-11}$ wherever M2 converges ($q\ge0.085$). Results are in METHODOLOGY §6.6a:

* thick-wall laws;
* error thresholds;
* finiteness at $13/72$ and $3/32$;
* crossovers;
* thin-wall comparison.

**SM lines.** $M_t=M_t^{\rm crit}+\delta$ with $\delta=10^{-4}$–$3$ GeV (17 points per line). $M_t^{\rm crit}$ is
where $\min\lambda_P=0$:

* 168.5030 GeV for (`cw1`, 120);
* 171.1193 GeV for (`cw1`, 125.15);
* 172.5079 GeV for (`cw1`, 128);
* 171.0689 GeV for (`tree`, 125.15).

| line | solved | $q/\delta$ [GeV$^{-1}$] | $\delta$ [GeV] where error = 1%: $S_{\rm FL}$ / $S_{\rm curv}$ / $S_{\rm CGM}$ / $S_{\rm CGM2}$ | error = 10%: same order | $S_{\rm CGM}$, $S_{\rm CGM2}$ worse than $S_{\rm FL}$ below $\delta$ = | unresolved |
|---|---|---|---|---|---|---|
| `cw1`, 120 | 15/17 | 121 | 0.152 / 0.0154 / 0.0175 / 0.0081 | 0.0132 / 0.0052 / 0.0054 / 0.0026 | 0.0034, 0.0014 | $\delta\le1.9\times10^{-4}$ ($q\le0.023$) |
| `cw1`, 125.15 | 14/17 | 141 | 0.130 / 0.0133 / 0.0149 / 0.0069 | 0.0113 / 0.0044 / 0.0046 / 0.0022 | 0.0029, 0.0013 | $\delta\le3.6\times10^{-4}$ ($q\le0.051$) |
| `tree`, 125.15 | 14/17 | 140 | 0.131 / 0.0134 / 0.0151 / 0.0070 | 0.0114 / 0.0045 / 0.0047 / 0.0023 | 0.0029, 0.0013 | $\delta\le3.6\times10^{-4}$ ($q\le0.051$) |
| `cw1`, 128 | 2/17 | – | – | – | – | 15 points with no finite bounce ($\delta\le0.83$ GeV) |

* **One parameter controls the regime.** Near criticality $q\propto\delta$, because $\lambda_{\min}\propto\delta$ at
  fixed curvature. On the solved points, $S/S_{\rm FL}-1$ follows the parabolic-model curve at the same $q$ to a
  median 2.2–2.5%. Over the controlled scans (standard tier, not Planck-sensitive; $q\in[0.17,600]$) the
  deviation is median 2.5% and at most 2.9%. The thresholds of the model therefore transfer to the SM through
  $q/\delta$. For example, the FL error reaches 1% at $q=17.8$ (model), which is $\delta=0.13$ GeV at
  $M_h=125.15$ GeV.
* **Scaling with $\delta$.** $d\ln S/d\ln\delta=-1.01$ at $\delta=1$–$3$ GeV, i.e. $S\propto1/|\lambda_{\min}|$
  (thick wall). It steepens to $-1.43$ (`cw1` and `tree`, 125.15) and $-1.69$ (`cw1`, 120) at the smallest
  solved $\delta$, as the wall thins.
* **The variational poles appear in the SM.** At $q=0.186$ (`cw1`, 125.15) $S_{\rm CGM}$ is 58 times too large
  and $S_{\rm curv}$ 16 times. At $q=0.184$ (`tree`) $S_{\rm CGM}=\infty$. $S_{\rm CGM2}$ is $\infty$ at $q=0.084$
  and 37–131 times too large at $q\approx0.097$, just above $3/32$.
* **Planck sensitivity.** Along `cw1`, $M_h=128$ GeV, the minimum of $\lambda$ lies at $10^{18.89}$–$10^{18.98}$
  GeV for $\delta\le0.83$ GeV, and no finite bounce exists (§5.2). The near-critical regime of that line is
  entirely Planck-dominated, so no statement about the analytical estimates can be made there. The other three
  lines are free of Planck sensitivity: their minima lie at $10^{17.1}$–$10^{18.4}$ GeV.
* **Unresolved points are genuine numerical failures.** For $q\lesssim0.05$, i.e. $\delta\lesssim4\times10^{-4}$
  GeV and $S\gtrsim5\times10^7$, neither the collocation solver (its Newton iteration does not converge to a bounce) nor
  the shooting fallback converges. This regime is not Planck-sensitive. Under the semiclassical assumptions of
  §5.2, the barrier bound $S_{\rm barrier}\ge S_{\rm FL}>1.1\times10^7\gg S_{\rm crit}$ makes them metastable.
  Their actions are not quoted.
* **Practical relevance.** Every point where an estimate errs by more than 1% has $S>3\times10^4$, in the
  scans as on the lines, where
  $\wp\sim e^{-S}$ is immaterial. At the central inputs of Buttazzo et al. ($M_t=173.34$ GeV, i.e.
  $\delta=2.2$ GeV at $M_h=125.15$ GeV) $q=258$, and $S_{\rm FL}$ is low by $7\times10^{-4}$.

### 5.3a All approximations on identical selections

Source: `validation/compare_approximations.py`, with output in `validation/results/approximation_comparison.json`
and `.md`. This is post-processing of the committed datasets only; no solver is run. It is tested in
`python/tests/test_approximation_comparison.py`, including a check that recomputation reproduces the committed
output.

* **Metric.** $e=(S_{\rm est}-S_{\rm num})/S_{\rm num}$. Medians and percentiles are empirical, without
  interpolation. An infinite or undefined estimate at a solved point counts as $|e|=\infty$.
* **Selection.** Within each dataset, every approximation is evaluated on the same points.
  * *Scans:* standard tier, i.e. $\lambda''>0$, not Planck-sensitive, not reduced-accuracy. For `cw1`, 10,432 of
    16,261 rows are included. The rest are 5,654 rows without a numerical action (stable, no finite bounce, or
    collocation failure) plus 174 Planck-sensitive and 1 reduced-accuracy point. For `tree`, 11,194 are included.
    For the exploratory scan, 1,699 of 10,100 are included.
  * *Near-critical SM lines:* all 45 solved points of 68. This includes 1 Planck-sensitive point (`cw1`, 128,
    $\delta=1.58$ GeV) and 4 reduced-accuracy points, among them the point at $q=0.044$ where $S_{\rm E1}$ has its
    largest error. Its numerical uncertainty, $\delta S/S=7\times10^{-6}$, is negligible against that error. The
    other 23 points are 15 without a finite bounce and 8 unresolved.
  * *Model family:* all 17 points.

  The full accounting is in the coverage table of `approximation_comparison.md`.
* **Approximations.**
  * $S_{\rm FL}$: rigorous lower bound.
  * $S_{\rm E1}=S_{\rm FL}(1+\frac{13}{72q})$: Espinosa (2020) eq. (3.21); a non-bound approximation.
  * $S_{\rm E2}$: adds Espinosa's second-order term ($c_2=-0.013644$); exact only for quadratic $\lambda$.
  * $S_{\rm curv}=S_{\rm FL}/(1-\frac{13}{72q})$: resummed first order; not a bound.
  * $S_{\rm CGM}$, $S_{\rm CGM2}$: variational upper estimates, conditional on CGM attainment.

**Controlled scans, both prescriptions (21,626 points).** Median $|e|$ / max $|e|$:

| region | points | $S_{\rm FL}$ | $S_{\rm E1}$ | $S_{\rm E2}$ | $S_{\rm curv}$ | $S_{\rm CGM}$ |
|---|---|---|---|---|---|---|
| $q\ge10$ (bulk) | 21,481 | $3.8\times10^{-4}$ / 0.018 | $9.2\times10^{-6}$ / $2.9\times10^{-4}$ | $9.3\times10^{-6}$ / $4.2\times10^{-4}$ | $9.1\times10^{-6}$ / $1.0\times10^{-4}$ | $2.1\times10^{-7}$ / $4.7\times10^{-4}$ |
| $1\le q<10$ | 118 | 0.025 / 0.144 | $3.0\times10^{-4}$ / $3.4\times10^{-3}$ | $5.7\times10^{-4}$ / $7.3\times10^{-3}$ | $3.4\times10^{-4}$ / 0.034 | $9.2\times10^{-4}$ / 0.039 |
| $0.1\le q<1$ (near-critical) | 27 | 0.226 / 0.504 | $8.5\times10^{-3}$ / **0.0154** | 0.021 / 0.212 | 0.110 / $\infty$ (1 point; largest finite 0.421) | 0.121 / $\infty$ (1 point; largest finite 0.459) |
| all | 21,626 | $3.9\times10^{-4}$ / 0.504 | $9.2\times10^{-6}$ / **0.0154** | $9.3\times10^{-6}$ / 0.212 | $9.1\times10^{-6}$ / $\infty$ | $2.1\times10^{-7}$ / $\infty$ |

Where the maxima occur:

* $S_{\rm E1}$: $+1.54\%$ at $(M_h,M_t)=(124.1,170.6)$ GeV, `cw1`, $q=0.338$, $S=2.6\times10^6$. The `tree` maximum is
  $+1.41\%$ at $q=0.173$.
* $S_{\rm FL}$: $-50.4\%$ at $(120.1,168.5)$ GeV, `tree`, $q=0.173$. $S_{\rm curv}$ and $S_{\rm CGM}$ are infinite at
  that point.
* The signed error of $S_{\rm E1}$ lies in $[-3.7\times10^{-4},+1.54\%]$.

**Other datasets** (same metric):

* *Exploratory scan* (1,699 standard-tier points, $q\ge2.9$): max $|e|$ is 6.0% ($S_{\rm FL}$), 0.23% ($S_{\rm E1}$),
  0.37% ($S_{\rm E2}$), 0.28% ($S_{\rm curv}$) and 0.57% ($S_{\rm CGM}$). The medians are $5.0\times10^{-4}$,
  $2.3\times10^{-5}$, $2.3\times10^{-5}$, $2.3\times10^{-5}$ and $3.6\times10^{-7}$.
* *SM near-critical lines* (45 solved points):
  * $0.1\le q<1$: max $|e|$ is 52% ($S_{\rm FL}$), 1.7% ($S_{\rm E1}$), 24% ($S_{\rm E2}$), $\infty$ ($S_{\rm curv}$,
    $S_{\rm CGM}$) and 56% ($S_{\rm CGM2}$).
  * $q<0.1$ (4 points): max $|e|$ is 83% ($S_{\rm FL}$), 14.9% ($S_{\rm E1}$, at $q=0.044$) and 133% ($S_{\rm E2}$).
    The variational estimates are infinite there, except $S_{\rm CGM2}$ at $q=0.097$ ($131\times$ too large).
* *Parabolic model* (17 points): for $0.1\le q<1$, max $|e|$ is 2.65% for $S_{\rm E1}$ (at $q=0.2$), 49% for
  $S_{\rm E2}$, $\infty$ for $S_{\rm curv}=S_{\rm CGM}$ and 712% for $S_{\rm CGM2}$. For $q<0.1$, $S_{\rm E1}$ reaches
  10.2% (at $q=0.05$).

**Findings.**

1. In the bulk, $S_{\rm CGM}$ has the smallest typical error, a median of $2\times10^{-7}$. That is about 40 times
   below the formulas that use only $\lambda''$ ($S_{\rm E1}$, $S_{\rm E2}$, $S_{\rm curv}$; about $10^{-5}$),
   whose residual comes from the cubic and higher terms of $\lambda(\ln\phi)$. $S_{\rm curv}$ has the smallest
   bulk maximum ($1.0\times10^{-4}$).
2. Near criticality, the unresummed first-order formula is by far the most accurate closed form: at most 1.54% on
   the scans, 1.7% on the SM lines for $q\ge0.097$, and 2.65% on the model for $q\ge0.1$. The resummed and
   variational estimates diverge there (METHODOLOGY §6.6a).
3. No estimate is uniformly best, and the accuracy of $S_{\rm E1}$ at $q\lesssim1$ is empirical, outside its
   formal domain. Below $q\approx0.07$ it degrades (−4% at $q=0.07$, −10% at 0.05 and −15% at the SM point
   $q=0.044$). The deep thin-wall regime $q\lesssim0.05$ is beyond both the formulas and the solvers.
4. These numbers come from three different domains that should not be merged: the SM scans ($q\ge0.17$), the SM
   lines ($q\ge0.044$) and the idealised parabolic model ($q\ge0.05$).

### 5.4 Exploratory wide scan (not validated)

`results/exploratory_cw1.csv` covers $[0,250]\times[2.5,250]$ GeV with extended matching fits, 10,100 points.
Its statuses:

| status | points |
|---|---|
| stable | 3817 |
| non-perturbative | 3545 |
| unstable | 1258 |
| no finite bounce | 811 |
| metastable | 538 |
| invalid input | 101 |
| solver failure | 30 |

815 points are flagged Planck-sensitive and 96 are in the reduced-accuracy tier. The solver failures were not
followed up. The matching fits are used far outside their domain, so the scan is an uncontrolled extrapolation
(`figures/exploratory_wide_scan.png`).

### 5.5 Phase-diagram scans (full plane and close-up)

The phase diagrams (README §10.7) use two further `cw1` scans with the extended matching: the full $[0,250]^2$ GeV
plane at 0.5 GeV (250,500 points) and the close-up $[110,140]\times[160,185]$ GeV at 0.05 GeV (301,101 points). The
scan outputs are not committed; their fallback results and `results/phase_diagrams_summary.json` are.

* **Agreement with the controlled scan.** At $M_h=125.2$ GeV the close-up gives $\lambda_{\min}=0$ at
  $M_t=171.144$ GeV and $\wp=1$ at 178.604 GeV. These are the controlled-scan values (README §10) to the quoted
  1 MeV.
* **Internal consistency.** The two scans give the same status at all 3,111 shared grid points, after the
  fallback and the action bound.
* **Recomputation.** A seeded, stratified sample, drawn as in `validation/reproduce_sample.py`, matched the stored
  values: 192 scan points recomputed from the clean build of `ad63002` at $10^{-12}$, and 16 fallback points
  re-solved (status, and $S$ at $10^{-10}$).
* **Solver failures.** Of 754 (full plane) and 4,433 (close-up) collocation failures, the shooting fallback decides
  730 and 4,417 with the acceptance criteria of §5.2. The action bound decides 23 and 7 more: a point with
  $\lambda_{\min}<0$ whose variational optimum lies at finite field is metastable if $S_{\rm FL}>559.3$, since then
  $\wp<1$ for every bounce scale up to $M_{\rm Pl}$. Ten points remain unresolved, all with the minimum of
  $\lambda$ at $10^{18.8}$–$10^{18.9}$ GeV.
* **Not validated.** Outside the controlled window the matching is extrapolated, as for the exploratory scan
  (§5.4).

## 6. Comparison with the legacy implementation

### 6.1 Why the legacy construction fails (precise mechanism)

The legacy `S_exact` inserts the Fubini–Lipatov profile with amplitude fixed by $\lambda_c\equiv\lambda(\mu=1/R)$
and integrates the potential over it. Algebraically this is identical to (METHODOLOGY §6.3)
$$S_{\rm legacy}(R)=\underbrace{\frac{8\pi^2}{3|\lambda_c|}}_{\propto 1/|\lambda_c|}
+\underbrace{\frac{\pi^2}{2}\int d\tau\,[\lambda(\phi_R)-\lambda_c]\,f_{\rm FL}^4}_{\propto\,\delta\lambda/\lambda_c^2},$$
which is the **first-order perturbative estimate** of the bounce action about the FL instanton of coupling
$\lambda_c$. The code then *minimises* it over $R$. It fails for four reasons:

1. **No validity check.** First-order perturbation theory requires $|\lambda-\lambda_c|\ll|\lambda_c|$ over
   the support of the profile. Nothing enforces this.
2. **Unbounded below.** The correction scales as $1/\lambda_c^2$ and the leading term as $1/|\lambda_c|$. Wherever
   $\lambda<\lambda_c$ somewhere in the profile while $\lambda_c\to0^-$, $S_{\rm legacy}\to-\infty$.
3. **Minimisation seeks out the invalid regime.** The bounce is a saddle point, not a minimum, so minimising over
   $R$ actively drives the construction into point 2. This happens with **perfectly smooth** potentials: on the
   parabolic model with $q=25$ the construction returns $-3\times10^{17}$, and near criticality in the SM
   $-1.3\times10^{21}$ (validation/results/method_comparison.csv; figures/legacy_failure.png panel a).
4. **Implementation defects that compound it** (SCIENTIFIC_METHODOLOGY F2, NUMERICAL_METHODS_REVIEW §3, §6):
   * an extra $\phi^6/(6M_{\rm Pl}^2)$ term;
   * $\lambda_c$ taken from the MS-bar coupling while the integrand uses $\lambda_{\rm eff}$;
   * piecewise-constant coupling interpolation, which makes the objective jagged at $10^{-6}$–$10^{-5}$;
   * a golden-section search with tolerance $10^{-13}$ on that jagged objective;
   * an early exit that returns $T$ instead of $T+V$, which created the artefact value 500000;
   * a $10^{100}$ sentinel that was classified as metastable.

On smooth potentials where the expansion *is* valid (thick wall), the construction is close to the true action but
still biased: by $10^{-4}$ at $q=250$, and by $-1.2\times10^{-4}$ at the SM central point. It is never a bound.

### 6.2 Legacy versus new on the controlled region

Source: `validation/results/legacy_comparison.csv` (`validation/legacy_comparison.py`). The legacy solver was
evaluated on a 0.5 GeV sub-grid (693 points) and compared with the new reference (`cw1`, Buttazzo matching):

| legacy \ new | metastable | unstable | stable | no-bounce | solver failure* |
|---|---|---|---|---|---|
| metastable | 276 | **56** | 0 | 7 | 15 |
| unstable | **8** | 98 | 0 | 0 | 0 |
| stable | **7** | 0 | 224 | 2 | 0 |

\* New-code collocation failures at Planck-sensitive points ($M_h\ge127$ GeV). After the shooting fallback (§5.2)
they are metastable or have no finite bounce (Planck plateau).

* **Agreement 86.3%.**
* **56 points: legacy metastable → new unstable.** The legacy threshold $4\ln(T\mu_{\rm inst})\approx485$ is 45–65
  units below the literature criterion.
* **8 points: legacy unstable → new metastable.** All eight have *negative* legacy actions.
* **7 points: legacy stable → new metastable.** These come from the corrected running: $g_2$ matching and the RGE
  transcription errors make $\lambda(M_{\rm Pl})$ more negative by $1.8\times10^{-3}$.
* **Where both are defined, the legacy $S_{\rm exact}$ exceeds the new reference by a median of +1.9%** (range −0.8%
  to +106%). This comes mostly from the $\phi^6$ term and the physics corrections.

The legacy outputs are preserved unchanged in `legacy/data/` and must not be used as physics results.

## 7. What is not validated

* **Physics beyond tree level.** No independent check of the decay-rate prefactor exists in this repository: the
  one-loop determinant is not computed.
* **Two-loop $\lambda_{\rm eff}$ and $e^{4\Gamma}$** are not implemented, so the absolute-stability boundary is not
  compared with Buttazzo et al. eq. (64) at their precision.
* **Inputs outside the controlled window** (the exploratory scan, and the phase diagrams beyond that window,
  §5.5) are not validated by construction.
* **Points flagged Planck-sensitive** depend on the arbitrary freezing of the couplings above $M_{\rm Pl}$. Their
  numerical solutions are accurate for that model, but the model itself is not physics. This includes the 459
  `cw1` points without a finite bounce.
* **The deep thin-wall regime** $q\lesssim0.05$ ($S\gtrsim5\times10^7$) is beyond the solvers (§5.3). There
  only the rigorous lower bound is available.
* **Solver agreement is not a proof.** The cross-method agreement of §3.2 is one line of evidence. The checks of
  §4 (self-convergence, precision, residuals, domain, quadrature, input perturbations) and the exact identities
  are independent of it. Even together they do not exclude an error common to the problem formulation, such as
  a wrong potential, which only the input checks of §1 address.

## 8. Uncertainty budget: numerical, physical and theoretical

These are kept separate. Only the first is an uncertainty *of this code*.

| kind | source | size | where |
|---|---|---|---|
| **numerical** | discretisation, round-off, domain, quadrature, input representation | $\sim10^{-14}$ relative (generic), $\sim10^{-12}$ (near criticality); per-point $\delta S\le1.1\times10^{-9}$ (standard tier); $\sim10^{-6}$ (reduced tier, Planck kink); shooting fallback $\le5.3\times10^{-10}$ | §4, §5.1–5.2 |
| **approximation** (closed-form vs numerical) | $S_{\rm FL}$, $S_{\rm E1}$, $S_{\rm E2}$, $S_{\rm curv}$, $S_{\rm CGM}$, $S_{\rm CGM2}$ | $3\times10^{-4}$–$2\times10^{-3}$ for $S_{\rm FL}$ in the bulk ($S<5000$); $\le1.54\%$ for $S_{\rm E1}$ over all standard-tier scan points; unbounded for the resummed and variational estimates near criticality | METHODOLOGY §6.5a, §6.6a; §5.3a |
| **parametric (physical inputs)** | $M_t$, $\alpha_s$, $M_h$ | not propagated here; the stability boundary moves by $\partial M_t^{\rm crit}/\partial M_h\approx0.5$ (from $M_t^{\rm crit}$ = 168.50, 171.12, 172.51 GeV at $M_h$ = 120, 125.15, 128) | §5.3 |
| **theoretical (RG improvement)** | prescription `cw1` vs `tree`; bounce scale $\Lambda_B=\arg\min\lambda$ vs $1/R$; omitted one-loop determinant, two-loop $\lambda_{\rm eff}$, $e^{4\Gamma}$ | 2.6–5% of $S$ (prescription); median 7.3 units in the decay threshold (scale); the determinant is not computed | §5.1, README §10.3 |
| **model (Planck scale)** | couplings frozen above $M_{\rm Pl}$; no gravity, no higher-dimensional operators | decisive for 459 `cw1` points (no finite bounce) and the 78 Planck-sensitive fallback bounces | §5.2 |
