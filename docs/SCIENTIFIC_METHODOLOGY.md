# Scientific Methodology — Reconstruction and Review

**Project:** SMVac — Standard Model electroweak vacuum stability
**Author:** Naman Goyal (GitHub: namangoyal31), Indian Institute of Science (IISc), Bangalore
**Status:** Stage B review document. It describes the physics implemented by the original
code (commit `a7c070a`, files `apps/solver_numerical.cpp` and the equivalent
`SMVacuumDecay/src/*`). Every statement is reconstructed from the source code, not from the
README. Where a statement depends on the literature, the source and equation are given.
Statements marked **[verified]** were checked numerically with the audit tooling in `audit/`.
Statements marked **[derived]** follow from a derivation given here.

Companion documents: `docs/REPOSITORY_AUDIT.md` (inventory, provenance, testing) and
`docs/NUMERICAL_METHODS_REVIEW.md` (algorithms, convergence, proposed replacements).

---

## 0. Summary of findings

| # | Finding | Category | Severity |
|---|---|---|---|
| F1 | The headline quantity `S_exact` is **not** a bounce action. It is a first-order perturbative estimate around the Fubini–Lipatov instanton, minimised over the instanton size outside its domain of validity. It produces negative "actions" at 8,740 grid points and the sentinel $10^{100}$ at 217 points, which are classified as *metastable*. | method | **fundamental** |
| F2 | A non-SM operator $V_6=\phi^6/(6M_{\rm Pl}^2)$ (coefficient 1, $M_{\rm Pl}=1.22\times10^{19}$ GeV) is added inside the `S_exact` integrand only. It is documented nowhere in this repository's README or code (only in the predecessor project's planning notes in `docs/`). At the SM point the whole 2.4% difference between `S_exact` and `S_approx` comes from this operator, not from RG effects. | physics | **fundamental for the headline claim** |
| F3 | The metastability threshold is evaluated at the *instability* scale (≈$10^{11}$ GeV) with $T\approx10$ Gyr. The literature criterion uses the *bounce* scale and the past-light-cone volume $0.15/H_0^4$. The legacy threshold is lower by 45–65 units, so the metastable region is too large. | physics | high |
| F4 | $g_2(M_t)=0.65355$ in the matching. The published NNLO value is 0.64779 (Buttazzo et al. 2013, eq. 58). This error alone shifts $\lambda(M_{\rm Pl})$ by $-1.5\times10^{-3}$, and the legacy value is 12% off. | physics input | high |
| F5 | Three transcription errors in the RG equations: the sign of the 4-loop QCD term, a missing $y_t^2$ in three 3-loop $\beta_\lambda$ terms, and a missing $y_t^2$ in a 3-loop $\beta_{y_t}$ group. | implementation | medium |
| F6 | The Z-boson term of $\lambda_{\rm eff}$ uses the GUT-normalised $g_1^2$ instead of $g_Y^2=\tfrac35 g_1^2$. This shifts $\lambda_{\rm eff,min}$ by $+4.4\times10^{-4}$, i.e. a 3.5% shift in the action. | implementation | medium |
| F7 | Linear matching fits valid near $(M_h,M_t)=(125,173)$ GeV are extrapolated over $[0,250]^2$ GeV, including $M_t\to0$ (negative $y_t$, $\ln 0$ starting scale). | validity | high (for the wide diagram) |
| F8 | $\lambda_{\rm eff}$ omits the anomalous-dimension factor $e^{4\Gamma}$, the Goldstone/Higgs terms and the two-loop correction that Buttazzo et al. include (their eqs. 104–107). | approximation | documented |

With F4–F6 corrected, the RG solution reproduces the published Planck-scale couplings of
Buttazzo et al. (2013), eq. (61), to the quoted precision **[verified]**. For the pure-SM
potential, a converged numerical solution of the tree-level bounce equation (computed here
independently) agrees with the simple Fubini–Lipatov estimate $8\pi^2/(3|\lambda_{\rm eff,min}|)$
to better than 0.1% at all benchmark points **[verified]**. That estimate is a rigorous lower
bound (§6.2).

**Fundamental problem (stop-and-explain item).** The repository's primary result is the phase
diagram with constant-`S_exact` contours, and its quoted validation numbers (the "32.1%
deviation") rest on F1 and F2. They are not supported by the implementation. The analytical
pipeline (`S_approx`) is a sound leading-order method once F3–F6 are corrected. No further
layers should be built on `S_exact`.

---

## 1. Physical problem and conventions

The Higgs doublet potential is $V=-\tfrac12 m^2|H|^2+\lambda|H|^4$, so at tree level
$M_h^2=2\lambda v^2$ with $v=246.22$ GeV. This is the normalisation of Buttazzo et al. (2013),
§3.1, and of the code's `lambda_tree = Mh^2/(2 v^2)`. For field values $h\gg v$ the potential is
$V(h)\simeq\lambda_{\rm eff}(h)\,h^4/4$. If $\lambda_{\rm eff}$ turns negative below the Planck
scale, the electroweak vacuum is not the global minimum. It is then *metastable* if its
lifetime exceeds the age of the Universe and *unstable* otherwise.

Conventions identified in the code **[verified against Buttazzo et al. App. B]**:

* RG time $t=\ln\bar\mu^2$; the code stores $t_0=2\ln M_t$.
* Gauge and Yukawa beta functions are $\beta_{x^2}=dx^2/d\ln\bar\mu^2$, with $LOOP_n=(16\pi^2)^n$.
  The integrator evolves $x$ via $dx/dt=\beta_{x^2}/(2x)$.
* $g_1$ is GUT-normalised, $g_1^2=\tfrac53 g_Y^2$ (one-loop coefficient $41/10$).
* MS-bar scheme throughout. The running starts at $\bar\mu=M_t$, with six active quark flavours.

---

## 2. Boundary conditions at $\bar\mu=M_t$ (`get_nnlo_matching`)

The code uses linear fits in $(M_h,M_t,\alpha_3)$ that mimic the NNLO results of Buttazzo et
al. (2013), §3, eqs. (55)–(60):

| coupling | code | Buttazzo et al. 2013 | status |
|---|---|---|---|
| $\lambda$ | $M_h^2/(2v^2)-0.00313-0.00004(M_t-173.34)$ | $0.12604+0.00206(M_h-125.15)-0.00004(M_t-173.34)$ (55) | equal at the central point (0.12605 vs 0.12604); tree-level quadratic $M_h$ dependence instead of a linear fit |
| $y_t$ | $\sqrt2 M_t/v-0.0587-0.00042\,\delta\alpha_3$ | $0.93690+0.00556(M_t-173.34)-0.00042\,\delta\alpha_3$ (57) | central value equal; slope 0.005744 vs 0.00556 (+3.3%) |
| $g_3$ | $1.1666+0.00314\,\delta\alpha_3-0.00046(M_t-173.34)$ | identical (60) | ✔ |
| $g_2$ | $0.65355+0.00002(M_t-173.34)$ | $0.64779+0.00004(M_t-173.34)+\dots$ (58) | **✘ central value wrong by +0.9%** (F4) |
| $g_1$ | $0.46266+0.00006(M_t-173.34)$ | $\sqrt{5/3}\,[0.35830+0.00011(M_t-173.34)]$ (59) | central value +0.02%, slope wrong by 2.4× |
| $y_b,\ y_\tau$ | $\sqrt2\,m/v$ with $m_b=4.0$, $m_\tau=1.777$ GeV | MS-bar values at $M_t$ (not given as fits) | tree-level with pole-like masses; impact negligible ($y_b^2\sim5\times10^{-4}$) |

Here $\delta\alpha_3=(\alpha_3(M_Z)-0.1184)/0.0007$. The code hard-codes $\alpha_3(M_Z)=0.1184$,
so $\delta\alpha_3\equiv0$. The fits are first-order expansions around
$(M_h,M_t)=(125.15,173.34)$ GeV. Their accuracy is controlled only within a few GeV of that
point. The original scan applies them over $M_h,M_t\in[0,250]$ GeV (F7):

* For $M_t<10.2$ GeV the fit gives $y_t<0$.
* For $M_t=0$ the starting scale is $\ln 0=-\infty$, and the integration loop terminates only
  when couplings blow up.
* For $M_t\lesssim M_Z$ the premise of matching to the full SM at $\bar\mu=M_t$, with six
  flavours, is invalid.

The wide-range diagram is therefore an extrapolation without theoretical control. Only the
neighbourhood of the measured masses is supported by these inputs.

## 3. Renormalisation-group evolution

**Implemented:** 3-loop gauge, $y_t$ and $\lambda$ beta functions; 2-loop $y_b$, $y_\tau$; the 4-loop
pure-QCD term in $\beta_{g_3^2}$. These are integrated from $M_t$ to $M_{\rm Pl}=1.22\times10^{19}$ GeV
(the non-reduced Planck mass), with no threshold below $M_t$ (none is needed above it).

**Coefficient verification [verified].** Every coefficient was compared term by term with
Buttazzo et al. (2013), Appendix B, eqs. (96)–(103), which use the same convention. All one-
and two-loop terms agree. Three discrepancies were found (F5):

1. **4-loop QCD sign.** The paper has $g_3^{10}(4\pi)^{-8}[-2472.28]$ (eq. 98); the code has
   $+2472.28$. The magnitude is the standard $\beta_3(n_f=6)$ of van Ritbergen, Vermaseren and
   Larin (1997), and in $d\alpha_s/d\ln\mu^2=-\sum\beta_i a^{i+2}$ form it enters with a minus sign.
2. **3-loop $\beta_\lambda$.** The paper has $y_t^6(-243.149y_t^2+250.494g_3^2+74.138g_2^2+33.930g_1^2)$
   (eq. 99). The code has $y_t^4(-243.149y_t^4+250.494g_3^2+\dots)$, so three terms lack a factor $y_t^2$.
3. **3-loop $\beta_{y_t^2}$.** The paper has $y_t^4(58.6028y_t^2+198\lambda-157g_3^2-\tfrac{1593}{16}g_2^2-\tfrac{2437}{80}g_1^2)$
   inside the bracket (eq. 101). The code has $y_t^2(\dots)$.

Items 2 and 3 are also detected by a literature-independent **loop-order power-counting test**
(`audit/check_power_counting.py`). Counting $g^2$, $y^2$ and $\lambda$ as degree one, every
$n$-loop term must have degree $n+1$. Exactly these two pieces violate it.

**End-to-end validation [verified]** (`audit/exp_defect_impact.py`; inputs
$M_h=125.15$, $M_t=173.34$, $\alpha_3=0.1184$; couplings at $1.22\times10^{19}$ GeV):

| | $g_1$ | $g_2$ | $g_3$ | $y_t$ | $\lambda$ |
|---|---|---|---|---|---|
| Buttazzo et al. eq. (61) | 0.6154 | 0.5055 | 0.4873 | 0.3825 | −0.0143 |
| legacy code | 0.6156 | 0.5082 | 0.4873 | 0.3810 | **−0.01257** |
| F4 + F5 corrected | 0.6154 | 0.5055 | 0.4873 | 0.3825 | −0.01437 |

Impact of each defect on its own (central point): $g_2$ matching
$\Delta\lambda(M_{\rm Pl})=-1.45\times10^{-3}$; 3-loop $\beta_\lambda$ $-2.7\times10^{-4}$; 3-loop
$\beta_{y_t}$ $-6\times10^{-5}$; 4-loop sign $-2\times10^{-5}$. Since
$\partial\lambda(M_{\rm Pl})/\partial M_t=-0.0066$ GeV$^{-1}$ (eq. 61), the combined error is
equivalent to a shift of about 0.27 GeV in $M_t$.

## 4. Effective potential (`get_lambda_eff`)

The code uses the one-loop Coleman–Weinberg correction in the large-field limit, evaluated at
$\bar\mu=h$ **[derived]**:
$$\lambda_{\rm eff}(h)=\lambda+\frac{1}{16\pi^2}\Big[-3y_t^4\big(\ln\tfrac{y_t^2}{2}-\tfrac32\big)
+\tfrac38g_2^4\big(\ln\tfrac{g_2^2}{4}-\tfrac56\big)+\tfrac{3}{16}(g_2^2+g_Y^2)^2\big(\ln\tfrac{g_2^2+g_Y^2}{4}-\tfrac56\big)\Big]_{\bar\mu=h}.$$
This follows from $V_1=\sum_i n_iM_i^4(\ln M_i^2/\bar\mu^2-c_i)/(64\pi^2)$ with
$n_t=-12$, $n_W=6$, $n_Z=3$, $c_t=3/2$ and $c_{W,Z}=5/6$ (Landau gauge, MS-bar). The top and W
terms match. **The Z term in the code uses $g_1^2+g_2^2$, where $g_1$ is GUT-normalised (F6).**
The correct combination is $g_Y^2+g_2^2=\tfrac35g_1^2+g_2^2$, as written in eq. (106) of Buttazzo
et al. Correcting it moves $\lambda_{\rm eff,min}$ from $-0.012903$ to $-0.012463$ at the
central point **[verified]**.

Pieces omitted relative to Buttazzo et al. eqs. (104)–(107) (F8; these are approximations to
document, not bugs):

* the field-renormalisation factor $e^{4\Gamma(h)}$ and the $2\Gamma$ shifts inside the logarithms;
* the Higgs/Goldstone term $3\lambda^2(4\ln\lambda-6+3\ln3)$, which is negligible near the instability;
* the two-loop correction $\lambda_{\rm eff}^{(2)}$;
* the quadratic term $-\tfrac12m^2h^2$, which is irrelevant for $h\gg v$.

The code sets $V=0$ below $h=246.22$ GeV, freezes all couplings above $M_{\rm Pl}$, and uses
the couplings at the lower grid node for $h$ between $v$ and $M_t$.

## 5. Stability classification

Status codes: 1 stable, 2 metastable, 3 unstable, 4 non-perturbative.

* **Stable:** $\lambda_{\rm eff}\ge0$ on every RG grid node from $\max(M_t,v)$ to $M_{\rm Pl}$.
* **Non-perturbative:** $|\lambda|>4\pi$ or $y_t>4\pi$ at any node. `FubiniLipatov.cpp` also uses
  $g_1>2$, so the two code paths use different criteria.
* **Metastable vs unstable:** decided by comparing an action with a threshold. `S_exact` is
  used in `classify_stability` and `S_approx` in `classify_buttazzo`; see §6–7.

The instability scale is recorded as the first grid node with $\lambda_{\rm eff}<0$ (Δln μ = 0.05
resolution). Using the non-reduced $M_{\rm Pl}$ as the upper end is a convention that should be
stated. Buttazzo et al. show results up to $M_{\rm Pl}\approx1.2\times10^{19}$ GeV.

## 6. Bounce action

### 6.1 Fubini–Lipatov instanton

For $V=\lambda\phi^4/4$ with constant $\lambda<0$, the O(4) bounce equation
$\phi''+3\phi'/r=V'(\phi)$ has the scale-invariant family
$\phi_R(r)=\sqrt{8/|\lambda|}\,R/(r^2+R^2)$ (Fubini 1976; Lipatov 1977). Its kinetic integral is
$T=16\pi^2/(3|\lambda|)$, its potential integral is $-8\pi^2/(3|\lambda|)$, and its action is
$S_0=8\pi^2/(3|\lambda|)$ for every $R$.

Two variables are useful throughout: $\tau=\ln r$ and $f=r\phi$. Then **[derived]**
$$S=2\pi^2\!\int d\tau\Big[\tfrac12f'^2+\tfrac12f^2+\tfrac14\lambda(\phi)f^4\Big],\qquad
f''=f+f^3\big(\lambda+\tfrac14\beta\big),\quad \beta\equiv\frac{d\lambda}{d\ln\phi},\quad \phi=fe^{-\tau}.$$
For constant $\lambda$ the solution is $f=\sqrt{2/|\lambda|}\,\mathrm{sech}(\tau-\tau_0)$; the translation
$\tau_0$ is the dilatation zero mode.

### 6.2 `S_approx` — the analytical estimate (sound)

`S_approx` $=8\pi^2/(3|\lambda_{\rm eff,min}|)$, with $\lambda_{\rm eff,min}$ the minimum over the RG
grid. This is the standard leading-order estimate.

**It is a rigorous lower bound** on the action of *every* finite-action O(4) bounce of
$V=\lambda(\phi)\phi^4/4$ **[derived]**. The proof:

* Any solution with finite kinetic energy satisfies the Pohozaev (Derrick) identity
  $T=-2V_{\rm int}$ in $d=4$, hence $S=T/2=T^2/(4|V_{\rm int}|)$.
* Pointwise $V\ge\lambda_{\min}\phi^4/4$ with $\lambda_{\min}<0$, so
  $|V_{\rm int}|\le\tfrac{|\lambda_{\min}|}{4}\int\phi^4$.
* The sharp Sobolev inequality in $\mathbb R^4$ (Aubin 1976; Talenti 1976) gives
  $\int\phi^4\le\frac{3}{32\pi^2}\big(\int|\nabla\phi|^2\big)^2=\frac{3}{8\pi^2}T^2$, with equality
  exactly for Fubini–Lipatov profiles.
* Combining: $S\ge8\pi^2/(3|\lambda_{\min}|)$.

Assumptions: $\lambda_{\min}$ is the global minimum of $\lambda(\phi)$ over all field values,
including any frozen regions, and the bounce has finite action. No existence or uniqueness
assumption is needed.

### 6.3 `S_exact` — what the original code actually computes (unsound)

`evaluate_action_at_R` inserts the FL profile $\phi_R$, with amplitude fixed by the MS-bar
coupling $\lambda_c\equiv\lambda(\bar\mu=1/R)$. It adds the exact kinetic term $16\pi^2/(3|\lambda_c|)$
and integrates $V=\lambda_{\rm eff}(\phi)\phi^4/4+\phi^6/(6M_{\rm Pl}^2)$ over the profile. Using
$\int f_{\rm FL}^4d\tau=16/(3\lambda_c^2)$, this is identically **[derived]**
$$S_{\rm legacy}(R)=\frac{8\pi^2}{3|\lambda_c|}+\frac{\pi^2}{2}\int d\tau\,\big[\tilde\lambda(\phi_R)-\lambda_c\big]f_{\rm FL}^4,
\qquad \tilde\lambda=\lambda_{\rm eff}+\tfrac23\phi^2/M_{\rm Pl}^2 .$$
This is the **first-order (Hellmann–Feynman) perturbative estimate** of the bounce action
around the FL instanton of coupling $\lambda_c$. The code then minimises it over $R$ by
golden-section search. Consequences:

* It is reliable only if $|\tilde\lambda-\lambda_c|\ll|\lambda_c|$ over the support of the profile.
  Nothing in the code checks this.
* The correction term scales like $\delta\lambda/\lambda_c^2$, while the leading term scales like
  $1/|\lambda_c|$. When $\lambda_c\to0^-$ at the edges of the region where the MS-bar $\lambda$
  is negative, any region with $\tilde\lambda<\lambda_c$ drives $S_{\rm legacy}\to-\infty$.
  **Minimisation over $R$ actively selects these invalid configurations.** This is the mechanism
  behind the 8,740 negative values, e.g. $S=-1.36\times10^6$ at $(M_h,M_t)=(5,105)$, and behind
  the "sandwiched unstable region" discussed in `docs/SANDWICH_*.md`.
* Where $\lambda_c\ge0$ the code returns the sentinel $10^{100}$. At 217 grid points the entire
  bracket returns it, and these points are classified *metastable*.
* Where $16\pi^2/(3|\lambda_c|)>5\times10^5$ the code returns the kinetic term alone ($=2S_0$).
  The search can converge onto this discontinuity: at $(5,96)$ the result is
  $500000.00000005$, and the plotting script adds a contour at exactly 500000.
* It is neither an upper nor a lower bound. At the SM point it exceeds the variational upper
  estimate for its own potential. At (115, 180) it falls below the rigorous lower bound (§6.6).

### 6.4 Lower bound and variational upper estimate (new, cheap) [derived, verified]

Coleman, Glaser and Martin (1978) characterise the least-action bounce in $d\ge3$ as the
minimiser of $T^2/(4|V_{\rm int}|)$ over configurations with $V_{\rm int}<0$. **The resulting upper
estimate is a rigorous bound only under the hypothesis that this infimum is attained by the
least-action bounce.** CGM establish attainment under conditions on $V$ that were not
re-verified for the RG-improved potential here, and attainment fails in the scale-invariant
degenerate case (§6.6, last row). We therefore call it a *variational upper estimate*. At every
benchmark point where a numerical bounce exists, we check numerically that it lies below this
estimate. Inserting the FL shape $\phi=h\,R^2/(r^2+R^2)$ and
changing variables to $w=R^2/(r^2+R^2)$ gives
$$\boxed{\;\frac{8\pi^2}{3\,|\min_\phi\lambda(\phi)|}\;\le\;S_b\;\le\;\frac{8\pi^2}{3\,|\min_h\bar\lambda(h)|},\qquad
\bar\lambda(h)=6\!\int_0^1\!w(1-w)\,\lambda(hw)\,dw\;}$$
The left inequality is rigorous (§6.2). The right one is conditional on CGM attainment. Both
sides are one-dimensional minimisations of smooth functions. They make no ODE solve and
no assumption of slow running. For constant $\lambda$ they coincide with $S_0$; their gap
measures how far the potential is from scale invariance over about one decade of field.

### 6.5 Numerical O(4) bounce (new) [verified]

The bounce of the RG-improved potential is a homoclinic orbit of the $f$ equation. The
mechanical energy $E=\tfrac12f'^2-\tfrac12f^2-\tfrac14\lambda f^4$ obeys
$dE/d\tau=\tfrac14\beta f^4$. Shooting on $\ln\phi(0)$ separates overshoot ($f$ crosses zero)
from undershoot ($f'$ returns to zero). A converged solution must satisfy two identities:

* the Derrick identity $T+2V_{\rm int}=0$;
* the dilatation identity $\int\beta f^4d\tau=0$, which is the Melnikov condition selecting $\tau_0$.

Implementation and tests: `audit/independent.py` and `audit/test_independent.py`. For constant
$\lambda$ the action is reproduced to $10^{-14}$. For a model with a smooth minimum of $\lambda$,
$S_{\rm lower}\le S_{\rm shoot}\le S_{\rm upper}$ holds, both identities hold to $10^{-13}$, and
the result converges under tolerance refinement.

### 6.6 Comparison at benchmark points [verified]

Same legacy physics inputs throughout (legacy matching, beta functions and $\lambda_{\rm eff}$),
so that only the bounce method differs (`audit/exp_bounce_comparison.py`):

| $(M_h,M_t)$ | potential | $S_{\rm lower}$ (= legacy `S_approx` for pure SM) | $S_{\rm shoot}$ (numerical) | $S_{\rm upper}$ (variational) | legacy `S_exact` |
|---|---|---|---|---|---|
| 125.1, 173.1 | pure SM | 2297.667 | **2299.595** | 2299.597 | 2353.818 |
| 125.1, 173.1 | SM + $\phi^6$ | 2343.93 | **2353.232** | 2353.285 | 2353.818 (> upper bound) |
| 120, 172 | pure SM | 1432.011 | **1432.881** | 1432.882 | 1435.609 |
| 120, 172 | SM + $\phi^6$ | 1434.01 | **1435.605** | 1435.608 | 1435.609 |
| 130, 178 | pure SM | 824.345 | **824.669** | 824.669 | 837.095 |
| 130, 178 | SM + $\phi^6$ | 835.38 | **837.209** | 837.214 | 837.095 |
| 115, 180 | pure SM ($\phi^6$ negligible here) | 327.463 | **327.5731** | 327.5731 | 327.4135 (**below the rigorous lower bound**) |
| 134.75, 176.5 | pure SM | 5740.921 | no finite bounce | 5740.921 | 7583.514 |

All shooting solutions satisfy the Derrick and dilatation identities to $\lesssim2\times10^{-12}$
(`audit/results/bounce_comparison.csv`). At (134.75, 176.5), $\lambda_{\rm eff}$ is still
decreasing at $M_{\rm Pl}$, where the couplings are frozen. The potential is therefore exactly
scale-invariant above $M_{\rm Pl}$. Every shot undershoots, and the two bounds coincide in an
infimum that no finite bounce attains: the minimising sequence escapes to arbitrarily large
field. The legacy 32% excess there is produced entirely by the $\phi^6$ term, which
regularises this Planck-scale degeneracy.
Conclusions:

* For the pure SM potential, the FL estimate is accurate to 0.04–0.08% and the upper bound is
  within $10^{-5}$ of the numerical action.
* The legacy `S_exact`−`S_approx` difference at these points is caused by the $\phi^6$ operator.
* Near the stability boundary, the README's 32% deviation at $(134.75,176.5)$ occurs where
  $\lambda_{\rm eff}$ is still decreasing at $M_{\rm Pl}$, the cut-off of the running. There the
  result is dominated by the Planck-scale treatment (couplings frozen, plus $\phi^6$), not by SM
  RG effects.

### 6.7 The hidden $\phi^6$ operator (F2)

`integrand_u` adds $V_6=\lambda_6\phi^6/(6M_{\rm Pl}^2)$ with $\lambda_6=1$. This is new physics
beyond the SM; Planck-suppressed operators are known to be able to change the vacuum lifetime
substantially (e.g. Branchina & Messina 2013). The operator appears only in `S_exact`. The
stability test and `S_approx` are pure SM, so the classification mixes two different theories.
If such operators are studied, they must be an explicit, documented model parameter
($\lambda_6$, $\Lambda$), and the cut-off treatment must be stated.

## 7. Decay probability and metastability criterion

Buttazzo et al. (2013), eqs. (80)–(82), give the decay probability in the past light cone:
$$d\wp=dt\,dV\,\Lambda_B^4e^{-S(\Lambda_B)},\qquad S(\Lambda_B)=\frac{8\pi^2}{3|\lambda(\Lambda_B)|},\qquad
\wp\simeq\frac{0.15}{H_0^4}\max_{\Lambda_B}\Lambda_B^4e^{-S(\Lambda_B)},\quad H_0\simeq67.4\ {\rm km\,s^{-1}Mpc^{-1}}.$$
Here $\Lambda_B=1/R$ is the bounce scale, in practice near the minimum of $\lambda$. The
boundary $\wp=1$ is therefore $S=4\ln(\Lambda_B/H_0)+\ln0.15$.

The code uses $S_{\rm thr}=4\ln(T\mu_{\rm inst})$ with $T=1.179\times10^{44}/v$ GeV$^{-1}$
$=9.99$ Gyr, not the cosmic age, and $\mu_{\rm inst}$ the instability scale (F3). At the SM point:

* legacy: $S_{\rm thr}=485.3$, with $\mu_{\rm inst}=1.0\times10^{11}$ GeV;
* literature: 530–550, with $\Lambda_B$ between $9\times10^{15}$ and $1.4\times10^{18}$ GeV.

The legacy threshold is too low by 45–65 units, which enlarges the metastable region. At the
measured masses $S\approx2000\gg550$, so the SM point is metastable under either criterion.

Two further limitations:

* The tree-level action ignores the one-loop fluctuation determinant and the related
  scale/scheme dependence. For example, using $\lambda_{\rm eff,min}$ instead of the MS-bar
  $\lambda_{\rm min}$ changes $S$ by about 10% at the SM point (`audit/results/defect_impact.csv`).
  This is a genuine theory systematic of any tree-level treatment.
* Gravitational corrections are ignored ($\delta S_G=256\pi^3\Lambda_B^2/(45|\lambda|M_{\rm Pl}^2)$,
  quoted in Buttazzo et al. §6.3).

## 8. Classification of assumptions

| Item | Type | Where |
|---|---|---|
| SM field content, MS-bar, Landau-gauge $\lambda_{\rm eff}$, no new physics below $M_{\rm Pl}$ | physical assumption | whole code |
| $V_6=\phi^6/(6M_{\rm Pl}^2)$ | physical assumption (BSM; not documented in README or code) | `integrand_u` |
| Large-field potential $\lambda_{\rm eff}h^4/4$; one-loop CW at $\bar\mu=h$; no $e^{4\Gamma}$ | analytical approximation | `get_lambda_eff` |
| Linear NNLO fits at $\bar\mu=M_t$ | analytical approximation (local) | `get_nnlo_matching` |
| Tree-level FL action; $\Gamma\propto\Lambda_B^4e^{-S}$ without prefactor | analytical approximation | `classify_*` |
| Couplings frozen above $M_{\rm Pl}$; $V=0$ below $v$ | modelling choice | `get_pure_sm_lambda` |
| RK4 with Δt = 0.1 in ln μ², grid-node minimum and zero crossing | numerical approximation (accurate, §NMR) | `classify_stability` |
| λ-only linear interpolation, other couplings piecewise constant | implementation choice (introduces roughness) | `RGEHelper` |
| Simpson rule with $N=2048$ on an algebraic map; golden-section search with tol $10^{-13}$, bracket ±8 in ln R, early exits | numerical method for an unsound objective | `evaluate_action_at_R`, `find_minimum_action` |

## References (verified to exist and to contain the cited content where marked)

1. D. Buttazzo, G. Degrassi, P.P. Giardino, G.F. Giudice, F. Sala, A. Salvio, A. Strumia,
   *Investigating the near-criticality of the Higgs boson*, JHEP 12 (2013) 089, arXiv:1307.3536.
   Eqs. (55)–(61), (80)–(82), (96)–(107) were read from the arXiv version for this review.
2. S. Fubini, Nuovo Cim. A 34 (1976) 521; L.N. Lipatov, Sov. Phys. JETP 45 (1977) 216.
3. S. Coleman, Phys. Rev. D 15 (1977) 2929; C.G. Callan, S. Coleman, Phys. Rev. D 16 (1977) 1762.
4. T. Aubin, J. Differential Geom. 11 (1976) 573; G. Talenti, Ann. Mat. Pura Appl. 110 (1976) 353
   (sharp Sobolev constant).
5. S. Coleman, V. Glaser, A. Martin, *Action minima among solutions to a class of Euclidean
   scalar field equations*, Commun. Math. Phys. 58 (1978) 211.
6. T. van Ritbergen, J.A.M. Vermaseren, S.A. Larin, Phys. Lett. B 400 (1997) 379,
   hep-ph/9701390. The value $\beta_3(n_f=6)=2472.28$ was recomputed from their formula.
7. G. Isidori, G. Ridolfi, A. Strumia, Nucl. Phys. B 609 (2001) 387, hep-ph/0104016.
   Cited for context; not re-read for this review.
8. V. Branchina, E. Messina, Phys. Rev. Lett. 111 (2013) 241801, arXiv:1307.5193.
   Cited for context; not re-read for this review.
9. A. Andreassen, W. Frost, M.D. Schwartz, Phys. Rev. D 97 (2018) 056006, arXiv:1707.08124.
   Cited for context; not re-read for this review.
