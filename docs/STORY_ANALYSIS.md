# STORY_ANALYSIS.md — Scientific Narrative Assessment

## Current Framing

The paper is currently framed as a **precision numerical methods paper** that happens to find a physical result. The title ("High-Precision Numerical Analysis...") emphasizes methodology over physics. The abstract spends the first 10 lines on numerics before stating the physics result.

## The Strongest Story

**"The SM vacuum stability boundary cannot be resolved by analytical approximations within current experimental uncertainties — exact numerical solvers are mandatory."**

This framing works because it:
1. **Creates tension** — the community consensus (conformal approximation) is insufficient
2. **Connects to experiment** — the 1σ M_t uncertainty matters
3. **Has a clear antagonist** — the Fubini-Lipatov / conformal approximation
4. **Offers resolution** — the exact numerical pipeline
5. **Generalizes** — the conclusion (exact solvers needed) applies beyond this specific calculation

## The Two Physical Discoveries

### Discovery 1: Barrier inflation rescues ~0.3–0.5 GeV of M_t
- Clean, intuitive physics: φ⁶ operator inflates the barrier → S_exact > S_approx
- Directly quantifiable: 110 parameter points reclassified
- Phenomenologically relevant: same size as 1σ(M_t)
- **This should be the paper's primary result.**

### Discovery 2: Fubini ansatz breakdown produces a "sandwiched" unstable region
- Deeper theoretical significance: reveals the regime of validity of the semi-classical approximation
- Virial ratio S_kin/|S_pot| ≠ 2 provides a quantitative diagnostic
- Explains an apparent paradox (negative action)
- Explains an apparent paradox (negative action) in a physically consistent way
- **This is the most interesting result for theory audiences but may confuse experimental readers.**
- **Recommendation:** Present it as a separate subsection with a clear pedagogical figure. If the paper targets PRD, keep it in the main text. If targeting PLB, relegate to an appendix.

## Narrative Arc (Recommended)

### Act 1: The Problem
- The SM vacuum is metastable (Buttazzo et al.)
- The conformal approximation is universally used but never tested near the boundary
- Near boundary, bounce profile reaches M_Pl → Planck operators matter

### Act 2: The Tool
- Build an exact numerical solver (summarize, don't dwell)
- Validate: four independent convergence tests → errors < 10⁻¹¹
- Scan 250,000 points in (M_t, M_h)

### Act 3: The Discovery
- Conformal approximation is fine globally (74.9% within 5%)...
- ...but fails catastrophically near the boundary (20–30% errors)
- φ⁶ operator rescues 0.3–0.5 GeV of parameter space
- Effect is universal across c₆

### Act 4: The Deeper Physics
- At low M_h, a sandwiched unstable region reveals Fubini ansatz breakdown
- Virial diagnostic identifies when the approximation fails
- Qualitatively new understanding of the semi-classical regime

### Act 5: The Implication
- Precision M_t measurements cannot resolve metastability vs. instability without exact solvers
- Roadmap for future work: include 1-loop CW corrections, temperature, UV completions

## Target Journals

| Journal | Fit | Concerns |
|---------|-----|----------|
| **Phys. Rev. D (PRD)** | Best fit — length, technical depth, topical | Length OK; sandwiched region may puzzle referees |
| **JHEP** | Good — similar scope to Buttazzo:2013uya | Length OK |
| **Phys. Lett. B (PLB)** | Too short for the technical content | Would need to cut most numerics |
| **Eur. Phys. J. C** | Good — phenomenological focus | Less prestige for pure SM vacuum stability |

## Comparison with Buttazzo:2013uya (The Canonical Reference)

| Aspect | Buttazzo et al. | This work |
|--------|-----------------|-----------|
| RG order | 3-loop | 3-loop |
| Bounce action | Conformal approx | Exact numerical |
| Planck ops | Not included | c₆ φ⁶/M_Pl² |
| Grid resolution | ~10² points | 2.5×10⁵ points |
| Error budget | Not quantified | < 10⁻¹¹ relative |
| Sandwiched region | Not seen | Identified and explained |

This comparison clarifies the niche: the paper is a **precision correction** to Buttazzo, not a replacement.

## Recommended Headline Numbers for the Abstract/Conclusion

1. **250,000** parameter points scanned
2. **110** points reclassified (unstable → metastable)
3. **~0.3–0.5 GeV** boundary shift in M_t (c₆ = 1)
4. **~1.5–2 GeV** shift for c₆ = 10
5. **< 10⁻¹¹** relative numerical error
6. **74.9%** of points within 5% of conformal approx (global accuracy)
7. **20–30%** discrepancy near boundary (failure regime)
