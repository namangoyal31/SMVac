# ABSTRACT_PLAN.md — Abstract Strategy

## Constraints

- **12–15 lines** at PRD two-column width (was 28 lines).
- No methodology details beyond one sentence.
- Lead with the physics problem, not the tool.
- State the headline number prominently.
- End with the implication.

---

## Three Versions

### Version A — "The Problem First" (Recommended)

Headline-first structure. Opens with the unresolved question, states the answer in line 2.

> The Standard Model electroweak vacuum may be metastable or unstable — a distinction that depends on the precise value of the top quark mass and on the accuracy of the semiclassical bounce action used to compute the tunneling rate. Standard analyses use the conformal approximation S ≃ 8π²/3|λ_min|, which ignores the shape of the potential at Planck-scale field values. By deploying an exact numerical solver (3-loop RGEs, RK4 integrator with ε = 10⁻¹³, Simpson quadrature N = 2048, golden-section radius optimization) across 250,000 (M_t, M_h) points, we show that this approximation systematically underestimates the bounce action by 20–30% near the metastability boundary. Including a Planck-scale φ⁶/M_Pl² operator with c₆ = 1 rescues a continuous strip of parameter space — approximately 110 points corresponding to 0.3–0.5 GeV in M_t — from unstable to metastable classification. The effect is universal for all c₆ > 0 and grows to 1.5–2 GeV for c₆ = 10. Since the boundary shift matches the 1σ experimental uncertainty on the top quark mass, exact numerical solvers are strictly necessary for precision vacuum stability determinations.

**12 lines.** Opens with tension, closes with implication.

### Version B — "The Surprise First"

Opens with the counterintuitive finding. Higher chance of grabbing a browser, higher risk of confusing.

> The conformal approximation for the electroweak vacuum decay rate — S ≃ 8π²/3|λ_min|, used in every major stability analysis since Buttazzo:2013uya — fails precisely where it matters most: near the metastability boundary, where it underestimates the true bounce action by 20–30%. An exact numerical computation across 250,000 (M_t, M_h) points reveals that including a Planck-scale φ⁶/M_Pl² operator (c₆ = 1) inflates the tunneling barrier enough to reclassify ∼110 parameter points (0.3–0.5 GeV in M_t) from unstable to metastable. The shift is comparable to the 1σ uncertainty on M_t ≈ 0.7 GeV, making exact solvers mandatory for precision vacuum stability analyses. The rescuing effect is universal for all c₆ > 0, reaching 1.5–2 GeV for c₆ = 10.

**10 lines.** Punchier, but skips the setup — non-specialist readers may not know what the conformal approximation is.

### Version C — "The Numbers First"

Least engaging, but maximizes the chance that a casual reader retains the headline number.

> We present an exact numerical computation of the electroweak vacuum decay bounce action across 250,000 points in the (M_t, M_h) plane. The three-loop RGE-improved potential includes a Planck-scale dimension-6 operator c₆ φ⁶/M_Pl². Total numerical errors are below 10⁻¹¹ in relative action. The conformal approximation S ≃ 8π²/3|λ_min| — accurate globally (74.9% within 5%) — underestimates the action by 20–30% near the metastability boundary. Including c₆ = 1 rescues ∼110 points (0.3–0.5 GeV in M_t) from unstable to metastable, a shift comparable to the 1σ uncertainty on the top quark mass (≈ 0.7 GeV). The rescuing is universal for all c₆ > 0, reaching 1.5–2 GeV for c₆ = 10. Exact numerical bounce solvers are therefore necessary for precision vacuum stability bounds.

**14 lines.** Slightly long, and the most important sentence ("the conformal approximation fails...") is buried. Not recommended.

---

## Element-by-Element Breakdown

| Element | Version A | Version B | Version C |
|---------|-----------|-----------|-----------|
| Opening hook | Problem stated as unresolved question | Counterintuitive finding | Our work / method |
| Method | 1 clause ("exact numerical solver...") | 1 clause | 2 clauses |
| Key number (0.3–0.5 GeV) | Line 8 | Line 5 | Line 8 |
| Implication | Last sentence | Last sentence | Last sentence |
| Tension level | High (question → answer) | Very high (surprise) | Low (method → result) |
| Risk | None | May confuse non-specialists | Buries the lead |

---

## Recommended Abstract Text

**Use Version A.** It respects the narrative arc (problem → tool → discovery → implication), works for both specialists and general readers, and puts the headline number (0.3–0.5 GeV) at exactly the right point in the paragraph.

---

## Title Alternatives

Ranked by estimated citation impact (1 = highest):

1. **Planck-Scale Rescuing of the Electroweak Vacuum** — concise, provocative, memorable.
2. **Exact Numerical Bounce Solver Reveals the Failure of the Conformal Approximation in SM Vacuum Stability** — descriptive, informative, but long.
3. **The Conformal Approximation is Insufficient: Exact Bounce Actions and Planck-Scale Operator Corrections to the SM Vacuum Stability Boundary** — precise but dry.
4. **Rescuing the Electroweak Vacuum: How Planck-Scale Operators Shift the Metastability Boundary by 0.5 GeV** — includes the headline number in the title (tactic used by high-impact papers).
5. **Beyond the Conformal Approximation: A 250,000-Point Exact Numerical Study of SM Vacuum Stability with Planck-Scale Operators** — emphasizes the scan size, signals thoroughness.

### My ranking

1. **Planck-Scale Rescuing of the Electroweak Vacuum** (option 1) — short, provocative, easy to cite and remember.
2. **Rescuing the Electroweak Vacuum: Planck-Scale Operators and the Failure of the Conformal Approximation** (hybrid of 1 + 2).
3. Option 4 (headline number in title).
4. Option 2 (descriptive).
5. Option 5 (numbers-focused).

**Recommendation:** Use a variant of option 1 or the hybrid. Avoid "High-Precision Numerical Analysis..." as a lead — it sounds like a methods paper.
