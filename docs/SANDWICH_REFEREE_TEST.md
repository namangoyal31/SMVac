# SANDWICH_REFEREE_TEST.md — Referee Verification Protocol for the Sandwiched Unstable Region

## Purpose

This document provides a skeptical referee with the minimal set of tests needed to independently verify the existence, characterization, and interpretation of the 2→3→2 sandwiched unstable region described in the paper.

---

## Test 1: Reproduce the Sandwich Pattern

### Minimal Verification
Run the following command to extract the 2→3→2 pattern at $M_h \approx 5$ GeV:

```python
import pandas as pd, numpy as np
df = pd.read_csv('results/numerical_data.csv')
mh_target, tol = 5.0, 0.02
chunk = df[(df['Mh_calc'] >= mh_target - tol) & 
           (df['Mh_calc'] <= mh_target + tol)].sort_values('Mt')
stabs = chunk['Stability'].values.astype(int)
mt_vals = chunk['Mt'].values
blocks = []
current_s = stabs[0]; current_start = mt_vals[0]
for i in range(1, len(stabs)):
    if stabs[i] != current_s:
        blocks.append((int(current_s), current_start, mt_vals[i-1]))
        current_s = stabs[i]; current_start = mt_vals[i]
blocks.append((int(current_s), current_start, mt_vals[-1]))
print('Pattern:', ''.join(str(b[0]) for b in blocks))
for b in blocks:
    print(f'  Stab={b[0]}: Mt [{b[1]:.1f}, {b[2]:.1f}] (width={b[2]-b[1]:.1f})')
```

**Expected output:** Pattern `412323234` with 2→3→2 at Mt ≈ [99, 114] GeV. The sandwiched-3 block should contain the Mt values [99.0, 99.5, ..., 113.8] with Stab=3.

### Test passes if:
- A 2→3→2 pattern is detected at M_h ≈ 5 GeV
- The sandwiched-3 block width is ~12–16 GeV
- Transition boundaries are sharp (no isolated single-point fluctuations)

---

## Test 2: Verify Negative-${S}$ in the Sandwiched Region

### Critical Test
The central claim of this analysis is that **100% of sandwiched-3 points have negative $S_{\rm exact}$**. Verify:

```python
# Continuing from Test 1
sw_mask = (chunk['Stability'] == 3) & (chunk['Mt'] >= 99) & (chunk['Mt'] <= 114)
sw_pts = chunk[sw_mask]
print(f'Points in sandwiched-3: {len(sw_pts)}')
print(f'Negative S_exact: {(sw_pts["S_exact"] < 0).sum()}')
print(f'Min S_exact in sandwich: {sw_pts["S_exact"].min():.1f}')
print(f'Max S_exact in sandwich: {sw_pts["S_exact"].max():.1f}')
```

**Expected output:**
```
Points in sandwiched-3: ~60
Negative S_exact: 60
Min S_exact in sandwich: < -1000000
Max S_exact in sandwich: < 0
```

### Test passes if:
- ALL points in the sandwiched-3 block have $S_{\rm exact} < 0$
- The negative values are large in magnitude ($10^3$–$10^7$), ruling out floating-point noise
- There is a discontinuity at the 3→2 (sandwich exit) boundary: e.g., from $S \approx -10^6$ directly to $S \approx +10^4$

---

## Test 3: Map the Full $M_h$ Extent of the Sandwich

Run `_sandwich_extent.py` to scan $M_h \in [0, 135]$ in 0.25 GeV steps:

```bash
python _sandwich_extent.py
```

**Expected output:** The 2→3→2 pattern should be present across the entire range $M_h \in [0, 135]$ GeV, though the sandwiched-3 width decreases from ~16 GeV at $M_h=0$ to ~0 GeV at $M_h\approx 100$, and becomes intermittent above 100 GeV.

### Test passes if:
- The sandwich exists at least for $M_h \lesssim 100$ GeV
- The total number of M_h slices with YES exceeds 200 (out of 540 scans)
- At M_h = 125 GeV (near the experimental point), there is NO sandwich
- The pattern chart shows clear 2→3→2 at low M_h and intermittent/absent at high M_h

---

## Test 4: Check for Negative Action Outside the Sandwich

Understand the scope of the Fubini-Lipatov breakdown:

```python
df = pd.read_csv('results/numerical_data.csv')
computed = df[df['S_exact'] != -1.0]
neg = computed[computed['S_exact'] < 0]
pos = computed[computed['S_exact'] > 0]
print(f'Computed points: {len(computed)}')
print(f'Negative S_exact: {len(neg)} ({100*len(neg)/len(computed):.2f}%)')
print(f'Min S_exact overall: {computed["S_exact"].min():.0f}')
print(f'Mt range of negative S: [{neg["Mt"].min():.1f}, {neg["Mt"].max():.1f}]')
print(f'Mh range of negative S: [{neg["Mh_calc"].min():.1f}, {neg["Mh_calc"].max():.1f}]')
```

### Test passes if:
- Negative $S$ is not confined to the sandwich — it appears broadly across the phase diagram
- This confirms the Fubini-Lipatov ansatz is unreliable wherever $S_{\rm exact}$ goes negative

---

## Test 5: Check for the Sandwich in the Dense Grid (if available)

If `results/numerical_dense_closeup.csv` exists, check for sandwich at $M_h \in [110, 140]$:

```python
d = pd.read_csv('results/numerical_dense_closeup.csv')
chunk = d[(d.iloc[:,1] >= 114) & (d.iloc[:,1] <= 116)].sort_values(d.columns[0])
stabs = chunk.iloc[:,2].values.astype(int)
# Check for 2→3→2 in the dense grid
unique_stabs = list(dict.fromkeys(stabs))
print(f'Stability sequence: {unique_stabs}')
```

**Expected result:** No 2→3→2 pattern at $M_h > 110$ GeV. The dense grid shows only a single 2→3 transition (the main metastability boundary).

### Test passes if:
- No 2→3→2 is found in the dense closeup region
- This confirms the sandwich is restricted to $M_h \lesssim 100$ GeV

---

## Test 6: Referee's Independent Consistency Check

The most definitive test would be to run the bounce solver with a **full ODE shooting method** (e.g., `SimpleBounce` or a custom implementation) at a sandwich point, e.g., $(M_t, M_h) = (105, 5)$ GeV where the Fubini-Lipatov ansatz gives $S_{\rm exact} < 0$.

### Predicted outcome:
- A full ODE solve would find a **positive** bounce action $S > 0$
- The value would likely be of order $S \sim \mathcal{O}(10^3)$ (metastable, not unstable)
- This would definitively confirm the sandwiched region is an ansatz artifact

This test requires external code (not in this repository) and is listed as a **recommended, not required, referee action**.

---

## Summary of Pass/Fail Criteria

| Test | What it proves | Pass condition |
|:----|:---------------|:--------------|
| 1 | Sandwich exists at M_h=5 | 2→3→2 pattern detected |
| 2 | All sandwich points have S<0 | 100% negative S |
| 3 | Full M_h extent mapped | Sandwich present across $\ge 200$ M_h slices |
| 4 | Negative S is broader than sandwich | $>15,000$ negative-S points total |
| 5 | Dense grid has no sandwich | No 2→3→2 at $M_h > 110$ |
| 6 | ODE solver would give S>0 | Predicted, not verifiable in repo |

---

## Verification Scripts in This Repository

| Script | Purpose |
|--------|---------|
| `_sandwich_map.py` | Rapid scan of all M_h for 2→3→2 pattern |
| `_sandwich_extent.py` | Detailed block-by-block analysis with S values |
| `_sandwich_detail.py` | S_exact profiling at specific M_h slices |
| `_sandwich_plot.py` | Generates all sandwich diagnostic figures |
| `_sandwich_check.py` | Quick M_h=5 check (from earlier audit) |

Run `python _sandwich_plot.py` to reproduce all figures in `results/`.
