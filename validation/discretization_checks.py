"""Discretisation, residual and round-off checks of the spectral reference solver (independent of solver agreement).

 1. Domain truncation: half-width L in {10, 12.5, 15, 20, 25} at fixed node density (N proportional to L).
 2. Off-node residual: the collocation residual vanishes at the nodes by construction, so the ODE residual
    f'' - f - f^3 G is evaluated on 20001 uniform points via the Chebyshev interpolant of the converged solution.
 3. Quadrature independence: the action recomputed with composite Simpson on the fine grid vs Clenshaw-Curtis.
 4. Round-off amplification: 10 runs with the Chebyshev coefficients of lambda(x) perturbed by random relative
    amounts of 2^-52 (one unit in the last place); the spread of S measures sensitivity to input round-off.
Output: validation/results/discretization_checks.txt
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from numpy.polynomial import chebyshev as C

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
from smvac.reference import analytic as A, bounce as B, potential as P, sm  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "results", "discretization_checks.txt")
lines = []


def say(s=""):
    print(s, flush=True)
    lines.append(s)


cases = [("parabola q=2.5", P.parabola_model(-2.5 * 5e-5, 5e-5, 39.0), (25, 55)),
         ("SM cw1 (125.15,173.34)", P.sm_potential(sm.Inputs(125.15, 173.34), "cw1"), None),
         ("SM cw1 (125.15,171.15)", P.sm_potential(sm.Inputs(125.15, 171.15), "cw1"), None)]

for name, pot, search in cases:
    an = A.estimates(pot, search)
    say(f"== {name}: S_FL={an.S_FL:.12g}")
    ref = None
    for L in (10.0, 12.5, 15.0, 20.0, 25.0):
        N = int(round(384 * L / 15.0 / 2)) * 2
        r = B.solve_m4(pot, an.lambda_min, an.lnh_CGM, N=N, L=L)
        if L == 15.0:
            ref = r
        say(f"   L={L:5.1f} N={N:4d}: S={r.S:.15g} derrick={r.derrick:.1e} dilatation={r.dilatation:.1e}")
    # off-node residual and alternative quadrature from the L = 15 solution
    tau, f = ref.profile
    tc, L = 0.5 * (tau[0] + tau[-1]), 0.5 * (tau[-1] - tau[0])
    xi = (tau - tc) / L
    coef = C.chebfit(xi, f, len(f) - 1)
    t = np.linspace(tau[0], tau[-1], 20001)
    u = (t - tc) / L
    F = C.chebval(u, coef)
    Fp = C.chebval(u, C.chebder(coef)) / L
    Fpp = C.chebval(u, C.chebder(coef, 2)) / L ** 2
    x = np.log(np.maximum(F, 1e-300)) - t
    G = pot.lam(x) + 0.25 * pot.dlam(x)
    res = Fpp - F - F ** 3 * G
    say(f"   off-node ODE residual: max|r|/max|f| = {np.max(np.abs(res)) / np.max(F):.2e}, "
        f"rms = {np.sqrt(np.mean(res ** 2)) / np.max(F):.2e}")
    h = t[1] - t[0]
    w = np.full(t.size, 2.0); w[1::2] = 4.0; w[0] = w[-1] = 1.0; w *= h / 3
    lam = pot.lam(x)
    T = 2 * math.pi ** 2 * (np.dot(w, 0.5 * (Fp - F) ** 2) + F[-1] ** 2)
    V = 2 * math.pi ** 2 * np.dot(w, 0.25 * lam * F ** 4)
    say(f"   Simpson (20001 pts) action: {T + V:.15g}  vs Clenshaw-Curtis {ref.S:.15g}  rel diff {abs((T + V) / ref.S - 1):.1e}")
    if hasattr(pot, "ser"):
        rng = np.random.default_rng(7)
        vals = []
        base = pot.ser.coef.copy()
        for k in range(10):
            pert = P.FrozenChebyshevPotential(lambda z: 0.0 * z, pot.x_lo, pot.x_hi, deg=8)
            pert.ser = C.Chebyshev(base * (1 + 2.0 ** -52 * rng.choice([-1, 1], size=base.size)), domain=[pot.x_lo, pot.x_hi])
            pert.d1, pert.d2 = pert.ser.deriv(1), pert.ser.deriv(2)
            pert.lam_lo, pert.lam_hi = float(pert.ser(pot.x_lo)), float(pert.ser(pot.x_hi))
            pert._tab = P.HermiteTable5(pot.x_lo, pot.x_hi, 5e-3, pert.ser, pert.d1, pert.d2)
            vals.append(B.solve_m4(pert, an.lambda_min, an.lnh_CGM, N=384).S)
        vals = np.array(vals)
        say(f"   ulp perturbation of lambda coefficients (10 runs): spread (max-min)/S = {(vals.max() - vals.min()) / vals.mean():.1e}")

with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("wrote", OUT)
