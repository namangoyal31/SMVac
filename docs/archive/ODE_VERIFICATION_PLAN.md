# ODE Verification Plan

**Test point:** M_h = 5 GeV, M_t = 105 GeV (deep inside the sandwiched-3 region)
**Diagnostic value at this point:** S_exact = -1.272x10^6, S_kin = 2.63x10^5, S_pot = -1.535x10^6, lambda_R = -2.00x10^-4, R_opt = 3.327, virial = 0.171

---

## 1. Which Bounce Equation Must Be Solved

The O(4)-symmetric Euclidean bounce action in flat space:

```
S[phi] = 2*pi^2 * integral_0^inf r^3 dr [ 0.5*(d phi/dr)^2 + V(phi) ]
```

The equation of motion (Euler-Lagrange) is:

```
d^2 phi/dr^2 + (3/r) d phi/dr = dV/d phi
```

with V(phi) given by the RG-improved SM effective potential using running couplings evaluated at mu = phi. The same SM RGEs as the existing solver are used (three-loop running for lambda, y_t, g1, g2, g3).

The equation is a second-order nonlinear ODE with a singular point at r = 0 (the 3/r term diverges). It must be integrated as an initial-value problem from r = 0 outward.

---

## 2. Boundary Conditions

### At r = 0 (center of the bounce)

```
phi'(0) = 0          (smoothness — O(4) symmetry)
phi(0)  = phi_0      (unknown — the shooting parameter)
```

Near r = 0, expand using L'Hopital's rule:

```
lim_{r->0} d^2 phi/dr^2 + (3/r) d phi/dr = 4 * d^2 phi/dr^2|_0

=>  d^2 phi/dr^2|_0 = (1/4) dV/d phi|_{phi_0}
```

Series expansion starting point:

```
phi(r) = phi_0 + (1/8) r^2 V'(phi_0) + O(r^4)
phi'(r) = (1/4) r V'(phi_0) + O(r^3)
```

Start numerical integration at r = r_0 = 10^-3 using the series — avoids the r=0 singularity.

### At r -> inf (far from the bounce)

```
phi(inf) = phi_FV ~ 0    (electroweak false vacuum)
phi'(inf) = 0
```

Asymptotic decay for a massive field:

```
phi(r) ~ A * exp(-sqrt(|V''(0)|) * r) / r^(3/2)    (for large r)
```

Terminate at r_max where |phi| < 10^-10.

### The shooting condition

Find phi_0 such that outward integration satisfies:

```
|phi(r_max)| < 10^-10    and    phi(r_max) > 0
```

---

## 3. Shooting Method Implementation

### Bisection loop

```
1. phi_lo = phi_FV + delta  (just above false vacuum — undershoot)
2. phi_hi = phi_TP          (turning point where V=0 — overshoot)
3. For each iteration:
   a. phi_0 = (phi_lo + phi_hi) / 2
   b. Integrate ODE from r=0 to r_max (RK4 adaptive)
   c. If phi crosses zero (overshoot): set phi_hi = phi_0
   d. If phi approaches false vacuum asymptotically: set phi_lo = phi_0
   e. If phi grows (overdamped): set phi_hi = phi_0
4. Repeat until |phi(r_max) - phi_FV| < 10^-10
```

Expected: 30-60 bisection iterations for double-precision accuracy.

### Integration technique

Rewrite as first-order system:

```
y1 = phi
y2 = d phi/dr

dy1/dr = y2
dy2/dr = -(3/r)*y2 + dV/d phi
```

Use RK4 with adaptive step size:
- h_min = 10^-4
- h_max = 1.0
- Local error tolerance: 10^-12
- Step size halving/doubling based on local truncation error estimate

### Termination logic

- phi < eps and dphi/dr < 0: converged (reached false vacuum)
- phi > phi_0 or dphi/dr > 0: undershoot (turning back)
- phi < -eps: overshoot (crossed zero)
- r > r_max: not converged (classify based on trend)

### Potential evaluation

Precompute V(phi) and dV/dphi(phi) on a log-spaced grid (10^2-10^9 GeV, 10^4 points) using the SM RGEs before integration begins. Use cubic spline interpolation during the ODE solve. This avoids re-running the RGE at every RK4 substep.

---

## 4. Numerical Difficulties Expected When lambda_R ~ 10^-4

### A. Extremely flat potential

When lambda_R ~ -2x10^-4, the quartic term V(phi) ~ lambda*phi^4/4 is nearly zero for moderate phi. The barrier is low and broad:

- dV/dphi is near zero over a large range of phi
- The restoring force is weak — the field rolls slowly
- Requires large r_max (potentially 10^4-10^7)
- The solution approaches the asymptotic form slowly

### B. Scale hierarchy may be manageable

The Fubini-Lipatov profile gives R_opt ~ 3.3 GeV^-1 and phi_0 ~ 100 GeV. If the true ODE bounce is near these scales, the integration is not pathological. However, if the true phi_0 is much larger (e.g., 10^6 GeV), a much wider radial domain is needed.

### C. Stiffness from RG running

lambda(mu) runs significantly across the integration domain (mu spans 10^2-10^6 GeV). Pre-tabulating V(phi) on a dense log grid and using spline interpolation handles this efficiently.

### D. Delicate shooting discriminant

When the potential is nearly flat, the difference between undershoot and overshoot is tiny. Mitigation:

- Tight ODE tolerance (10^-12)
- Verify convergence by doubling r_max and checking S changes < 10^-6
- Use high-precision (long double) for the ODE integration

### E. V_min sensitivity to M_t

At M_h = 5 GeV, the true vacuum depth depends on RG running at 10^6-10^10 GeV. Vary M_t by +/-0.7 GeV and compute S for each case to establish the physical uncertainty band.

---

## 5. How to Compare Against the Fubini-Lipatov Ansatz

### Direct S_total comparison

The FL ansatz gives S = -1.27x10^6 (negative, unphysical). The ODE should give S > 0.

```
Delta_S = S_ODE - S_FL
```

Classification implications:
- S_ODE > 400: point is metastable (FL false positive)
- 0 < S_ODE < 400: point is genuinely unstable (FL got sign wrong but classification correct by coincidence)
- S_ODE < 0: something wrong with ODE solver

### Profile shape comparison

FL profile:

```
phi_FL(r) = sqrt(2/|lambda_R|) * 2R / (r^2 + R^2)
```

Compare phi_ODE(r) vs phi_FL(r):
- Log-linear overlay plot
- L2 norm of the difference
- Core amplitude phi_ODE(0) vs sqrt(2/|lambda_R|)
- Asymptotic tail exponent

### Characteristic radius

Use R_half: phi(R_half) = phi(0)/2. Compare R_ODE vs R_FL = 3.33 GeV^-1. If they differ by >30%, the FL profile shape is wrong even with the optimal radius.

### Virial ratio

The most diagnostic comparison:

```
S_kin = 2*pi^2 * integral r^3 (0.5*(d phi/dr)^2) dr
S_pot = 2*pi^2 * integral r^3 V(phi) dr
```

FL at this point: S_kin/|S_pot| = 0.171
ODE expected: S_kin/|S_pot| ~ 2 (Derrick's theorem for exact O(4) bounce)

If the ODE virial ratio is ~2, this is definitive proof that the ansatz shape constraint (not the physics) causes the negative S.

---

## 6. Which Observables Should Be Compared

### Primary (must-report)

| Observable | FL value | ODE expected | What it tells you |
|-----------|:--------:|:------------:|-------------------|
| S_total | -1.27e6 | > 0 | Sign and magnitude |
| S_kin | 2.63e5 | TBD | Kinetic term |
| S_pot | -1.54e6 | TBD | Potential term |
| S_kin/|S_pot| | 0.171 | ~2 | Virial theorem: FL validity diagnostic |
| phi_0 = phi(0) | sqrt(2/|lambda_R|) ~ 100 GeV | TBD | Core field value |
| R_half | 3.33 GeV^-1 | TBD | Characteristic radius |

### Secondary (diagnostic)

| Observable | Purpose |
|-----------|---------|
| phi(r) profile | Visual shape comparison |
| d phi/dr (r) | Where kinetic energy is concentrated |
| phi(r) - phi_FL(r) | How far the ansatz deviates |
| S(total) vs integration limit | Convergence diagnostic |
| M_t +/- 0.7 GeV | Physical uncertainty on S |

### Classification outcome

| ODE S_total | Meaning |
|:-----------:|---------|
| > 400 | Metastable. FL was wrong — false positive for instability |
| 0-400 | Unstable. FL got the classification right despite negative S |
| < 0 | ODE solver error — unexpected |

---

## Summary: One-Point Validation Test

```
Point:       M_h = 5 GeV, M_t = 105 GeV
FL result:   S = -1.27e6, virial = 0.171, R = 3.33 GeV^-1

Procedure:
1. Precompute V(phi), dV/dphi on log-spaced grid (10^2-10^9 GeV, 10^4 pts)
2. Implement RK4 adaptive (h_min=1e-4, h_max=1.0, tol=1e-12)
3. Start at r=1e-3 using series expansion
4. Shoot phi_0 in [10, 10^4] GeV via bisection (50 iterations)
5. Terminate at r_max=10^4 or |phi|<1e-10
6. Compute S_kin, S_pot, virial, R_half, phi_0

Expected outcome:
- S_ODE > 0 (probably 10^2-10^4)
- virial ratio ~ 2 (confirming Derrick's theorem)
- R_ODE != R_FL (ansatz gets the scale wrong)
- Proves the sandwich is an FL ansatz artifact
```
