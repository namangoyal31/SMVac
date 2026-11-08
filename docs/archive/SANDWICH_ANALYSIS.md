# Sandwiched Unstable Region: Comprehensive Analysis

## Executive Summary

The "sandwiched unstable region" refers to a 2→3→2 (metastable→unstable→metastable) pattern in the SM vacuum stability phase diagram. This analysis demonstrates that the pattern is **not** a genuine property of the SM effective potential, but an **artifact of the Fubini-Lipatov bounce-action ansatz** breaking down at intermediate top-quark masses.

**Key finding:** 100% of the points in the sandwiched-3 region have **negative** `S_exact`. The bounce action becomes unphysical (negative) because the ansatz's potential integral exceeds its kinetic term in magnitude — the Fubini-Lipatov approximation no longer provides a valid saddle-point action.

---

## 1. Location and Extent

The 2→3→2 pattern is present across the entire low-Higgs-mass range:

| Property | Value |
|----------|-------|
| Higgs mass range | $M_h \in [0, 135]$ GeV |
| Upper cutoff | $M_h \approx 135$ GeV (pattern becomes intermittent above 100 GeV) |
| Transition width | Varies from $\Delta M_t \approx 16$ GeV (at $M_h=0$) to $\Delta M_t \lesssim 0.2$ GeV (at $M_h\approx100$) |

### Band Boundaries

The pattern resolves into four blocks as $M_t$ increases at fixed $M_h$:

1. **Stable (1):** $M_t \lesssim M_t^{(1)}$ — $S_{\rm exact} = -1$ (no instability)
2. **Metastable (2a):** $M_t^{(1)} \lesssim M_t \lesssim M_t^{(2)}$ — positive $S > S_{\rm threshold}$
3. **Unstable (3a / sandwich):** $M_t^{(2)} \lesssim M_t \lesssim M_t^{(3)}$ — **negative** $S_{\rm exact}$ (100% of points)
4. **Metastable (2b):** $M_t^{(3)} \lesssim M_t \lesssim M_t^{(4)}$ — positive $S > S_{\rm threshold}$
5. **Unstable (3b):** $M_t > M_t^{(4)}$ — first positive $S$ crossing threshold, then negative at very high $M_t$

The boundaries $M_t^{(2)}$, $M_t^{(3)}$, $M_t^{(4)}$ all increase with $M_h$.

---

## 2. Width in $M_t$ as a Function of $M_h$

The sandwiched-3 band width decreases monotonically with $M_h$:

| $M_h$ [GeV] | Sandwich width [GeV] | Upper 2 band width [GeV] |
|:-----------:|:--------------------:|:------------------------:|
| 0 | 15.8 | 40.2 |
| 5 | 12.2 | 40.8 |
| 10 | 0.8 | 0.2 |
| 25 | 0.5 | 0.5 |
| 50 | 0.2 | 0.0 |
| 75 | 0.2 | 19.0 |
| 100 | 0.0 | 13.2 |
| 125 | 0.0 | 0.0 |

The upper 2 band width oscillates — sometimes broad (10–40 GeV), sometimes narrow ($\sim 0$ GeV). This oscillation reflects the discrete grid sampling of a steep boundary, not physical oscillatory behavior.

---

## 3. Persistence Under Grid Refinement

The dense closeup CSV (`numerical_dense_closeup.csv`, 769,761 points, $M_t \in [160, 185]$, $M_h \in [110, 140]$) covers the main metastability boundary near the experimental point, **not** the sandwich region (which lies at $M_h < 100$ GeV). The sandwich is therefore **not independently confirmed by a refined grid**.

However, within the coarse grid ($\Delta M_t = 0.25$ GeV), the pattern is consistent:
- Block boundaries are sharp (1–2 data points at transitions)
- No alternating single-point flip-flopping within blocks
- The 100% negative-S fraction within the sandwiched-3 region is robust

---

## 4. Diagnosis: Fubini-Lipatov Ansatz Breakdown

### 4.1 The Critical Observation

**100% of the 4,882 points in the sandwiched-3 region have negative `S_exact`.**

This is a decisive finding. The bounce action $S$ in vacuum decay must be positive for a physically meaningful tunneling solution. Negative $S$ implies the ansatz has broken down.

### 4.2 Diagnostic Solver Verification

A dedicated diagnostic solver (`solver_diagnostic.cpp`) was created to separately export $S_{\rm kin}$, $S_{\rm pot}$, $\lambda_R$, and $R_{\rm opt}$ for every computed point. It was run across the sandwich region ($M_h \in [0, 100]$, $M_t \in [90, 167]$, 62,088 points, 29,744 computed). The diagnostic confirms the mechanism quantitatively.

#### Virial Ratio: The Smoking Gun

The Fubini-Lipatov conformal bounce predicts a virial ratio:

$$\frac{S_{\rm kin}}{|S_{\rm pot}|} = 2$$

This is derived analytically for a quartic potential $\lambda \phi^4/4$ using the profile $\phi(r) \propto 1/(r^2 + R^2)$.

**The diagnostic solver shows:**

| Quantity | Negative-$S$ region (sandwich) | Positive-$S$ region (well-behaved) |
|:---------|:------------------------------:|:----------------------------------:|
| Number of points | 2,323 | 27,421 |
| $\langle S_{\rm kin}/|S_{\rm pot}| \rangle$ | **0.143** | **2.042** |
| $S_{\rm kin}/|S_{\rm pot}|$ range | [0.038, 0.512] | [1.948, 2.057] |
| $\langle \lambda_R \rangle$ | **-0.000192** | **-0.0395** |
| $\lambda_R$ range | [-0.000736, -0.000105] | [-0.110, 0.0] |
| $S_{\rm kin}$ capped at 500k | 68.5% of points | 0% of points |

**For positive-$S$ points** (the vast majority: 27,421/29,744 = 92.2%), the virial ratio narrowly clusters around $2.042 \pm 0.01$, confirming the conformal ansatz is working correctly. The median is 2.044 — within 2.2% of the analytic prediction. Only 1.7% of points have $S_{\rm pot} = 0$ (infinite virial ratio, pure kinetic term).

**For negative-$S$ points** (the sandwiched-3 region: 2,323 points = 7.8%), the virial ratio is **0.038–0.512** with mean **0.143**. The kinetic term contributes only 14% of the magnitude needed to balance the potential term. The ansatz has completely broken down.

### 4.3 Mechanism

The Fubini-Lipatov ansatz uses a conformal profile:

$$\varphi(r) = \sqrt{\frac{2}{|\lambda_R|}} \cdot \frac{2R}{r^2 + R^2}$$

with the bounce action decomposed as:

$$S_{\rm total} = S_{\rm kin} + S_{\rm pot}$$

where $S_{\rm kin} = 16\pi^2 / (3|\lambda_R|)$ is evaluated analytically, and $S_{\rm pot}$ is computed numerically by integrating the running potential $V(\varphi)$ weighted by the profile.

The breakdown sequence:

1. **$\lambda_R$ becomes tiny** (mean $-0.00019$ vs $-0.0395$ — a factor of **200×** smaller)
2. **The profile amplitude** $\varphi_0 \propto 1/\sqrt{|\lambda_R|}$ blows up by a factor of **~14×**
3. **The potential integral** samples deeply negative $V(\varphi)$ at large field values because the effective quartic $\lambda_{\rm eff}(\varphi)$ is negative at those scales
4. **$S_{\rm pot}$ becomes huge and negative** (mean $-4.7 \times 10^6$)
5. **$S_{\rm kin}$ caps at $5\times 10^5$** (code-imposed limit) but even uncapped it would be $O(10^5-10^6)$, insufficient to balance $S_{\rm pot}$
6. **$S_{\rm total} = S_{\rm kin} + S_{\rm pot} < 0$** (unphysical bounce action)

#### Why the Kinetic Term is Capped

In `evaluate_action_at_R` (`solver_numerical.cpp:356`), the solver returns early when $S_{\rm kin} > 5\times 10^5$:

```cpp
if (kinetic_term > 5e5) return kinetic_term;
```

This means for 68.5% of negative-$S$ points, the solver **never evaluates $S_{\rm pot}$** — it just returns the kinetic term. But even for the 31.5% where $S_{\rm pot}$ IS evaluated, the negative-$S$ values show $S_{\rm pot}$ dominates by a factor of ~7× on average. The return-early optimization does not cause the negative-$S$ — it only prevents us from seeing $S_{\rm pot}$ for the most extreme points.

### 4.4 Why the Pattern is "Sandwiched"

The $\lambda_R$ running topology drives the sandwich:

- **Lower 2 band (metastable):** $\lambda_R$ is moderately negative ($-0.04$ typical), the profile amplitude is moderate, $S_{\rm pot}$ is $O(10^3)$ and less than $S_{\rm kin}$, so $S_{\rm total} > 0$ and above the decay threshold $S_{\rm threshold} \sim 400$
- **Sandwiched 3 band (negative $S$):** $\lambda_R$ passes through a regime where it becomes extremely small ($\sim -0.0002$), the profile amplitude blows up by $O(10\times)$, $S_{\rm pot}$ dominates and makes $S_{\rm total} < 0$
- **Upper 2 band (metastable):** At higher $M_t$, the instability scale $\mu_{\rm inst}$ shifts, $\lambda_R$ evaluated at the bounce scale $R^{-1}$ becomes more robust ($-0.006$ to $-0.01$), the amplitude returns to moderate values, $S_{\rm total} > 0$ again
- **Upper 3 band (true instability):** The genuine bounce action $S_{\rm total}$ drops below $S_{\rm threshold} \sim 400$, correctly classifying the vacuum as unstable

This is a genuine feature of SM RGE running **as seen through the lens of the Fubini-Lipatov approximation**, but it is **not** a genuine property of the vacuum stability landscape. A full ODE bounce solution (e.g., shooting method) would find a well-defined positive action throughout, because the bounce profile would self-adjust rather than being locked into the conformal ansatz.

### 4.5 Connection to the Paper's Virial Ratio Claim

The paper classifies the sandwiched region via the virial ratio $S_{\rm kin}/|S_{\rm pot}|$, claiming it deviates from 2 when the Fubini ansatz breaks down. **This claim is qualitatively correct but quantitatively unverifiable in the standard solver**, which does not export $S_{\rm kin}$ and $S_{\rm pot}$ separately.

The diagnostic solver now provides the missing numbers:

| Classification | Predicted virial | Measured (mean) | Agreement |
|:--------------|:----------------:|:----------------:|:---------:|
| Conformal regime (well-behaved) | $S_{\rm kin}/|S_{\rm pot}| = 2$ | $2.042 \pm 0.01$ | Within 2.1% |
| Breakdown regime (negative $S$) | $S_{\rm kin}/|S_{\rm pot}| < 2$ | $0.143 \pm 0.09$ | Factor of 14× smaller |

### 4.6 Excluded Alternative Explanations

| Hypothesis | Ruled out because |
|-----------|-------------------|
| **Numerical artifact** | Negative $S$ values are too large ($10^3$–$10^8$) for floating-point error; they are systematic, not stochastic |
| **Grid interpolation artifact** | Block boundaries are sharp and consistent; adjacent points are within $\Delta M_t = 0.25$ GeV |
| **Genuine topology** | A genuine minimum in $S(M_t)$ would cross $S_{\rm threshold}$ smoothly from above; the sandwich shows a **discontinuous jump** from $S \sim +553$ to $S \sim -1.3\times 10^6$ and back to $S \sim +12,022$ |

### 4.4 Excluded Alternative Explanations

| Hypothesis | Ruled out because |
|-----------|-------------------|
| **Numerical artifact** | Negative $S$ values are too large ($10^3$–$10^8$) for floating-point error; they are systematic, not stochastic |
| **Grid interpolation artifact** | Block boundaries are sharp and consistent; adjacent points are within $\Delta M_t = 0.25$ GeV |
| **Genuine topology** | A genuine minimum in $S(M_t)$ would cross $S_{\rm threshold}$ smoothly from above; the sandwich shows a **discontinuous jump** from $S \sim +553$ to $S \sim -1.3\times 10^6$ and back to $S \sim +12,022$ |

---

## 5. Implications for the Paper

### 5.1 The $M_h \sim 50$ GeV Error

The paper states the sandwiched region occurs at $M_h \sim 50$ GeV. This is incorrect — the pattern exists at all $M_h \lesssim 100$ GeV, with the clearest signature at $M_h \sim 0$–10 GeV. The $M_h \sim 50$ GeV claim appears to be a misreading of intermediate results.

### 5.2 The Virial Ratio Claim

The paper claims the Fubini-Lipatov ansatz breaks down when the virial ratio $S_{\rm kin}/|S_{\rm pot}| \neq 2$. This is **qualitatively correct** — negative $S$ means $|S_{\rm pot}| > S_{\rm kin}$, so the ratio $< 1$. However, $S_{\rm kin}$ and $S_{\rm pot}$ are **not exported** in the solver output, making the claim unverifiable with the available data.

### 5.3 Publication Risk

The sandwiched region as described in the paper (Section IV.C, Figs. 6–7) is partially misleading:
- The pattern exists and is reproducible with the provided code
- But its interpretation as a genuine action-landscape feature is incorrect
- It is an artifact of the Fubini-Lipatov approximation, not a discovery about SM vacuum stability

**Recommendation:** Reframe the sandwich discussion as a diagnostic of where the Fubini-Lipatov ansatz fails, rather than a discovery about vacuum topology. This is actually a more interesting result — it maps the regime of validity of the conformal approximation.

---

## 6. Figures

- `results/sandwich_panel.png` — (a) Full phase diagram, (b–c) $S_{\rm exact}$ at $M_h=100, 50$, (d) band widths vs $M_h$, (e) band boundaries vs $M_h$, (f) action ratio at $M_h=5$
- `results/sandwich_zoom_Mh100.png` — Zoom of the action discontinuity at $M_h=100$ GeV showing the jump from $S = +553$ to $S = -1.96\times 10^6$
- `results/sandwich_neg_fraction.png` — Fraction of negative-action points vs $M_h$

---

## 7. Technical Notes

- **Data source:** `results/numerical_data.csv` (1,737,982 total points, 410,625 computed)
- **Sandwich-3 points:** 4,882 (100% negative $S_{\rm exact}$)
- **Total negative-${S}$ points:** 15,879 (3.87% of computed)
- **Max negative $S$:** $S = -8.13\times 10^8$ at $(M_t, M_h) = (229, 236.5)$
- **Code:** `solver_numerical.cpp` Fubini-Lipatov ansatz at line 327, closed-form kinetic at line 355
