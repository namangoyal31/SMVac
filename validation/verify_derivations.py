"""Independent numerical checks of the analytical results used in SMVac (docs/SCIENTIFIC_METHODOLOGY.md section 6).

These checks do not prove the theorems; they test the constants, formulas and stated orderings against direct
computations that share no code with the production implementation.

 1. Sharp Sobolev constant in R^4: K = 3/(32 pi^2) is attained by the Fubini-Lipatov (FL) profile and not exceeded
    by other radial trial profiles  -> the ingredient of the rigorous lower bound S >= 8 pi^2/(3|lambda_min|).
 2. Variational formula: T^2/(4|V|) evaluated by direct quadrature in r for the trial phi_h(r) = h/(1+r^2) equals
    (8 pi^2/3)/|lambdabar(h)| with lambdabar(h) = 6 int_0^1 w(1-w) lambda(h w) dw.
 3. Curvature coefficient: E[ln w] = -5/6 and (1/2) Var[ln w] = 13/72 under the FL weight 6 w (1-w) dw.
 4. First-order asymptotics: for lambda = l0 + (l2/2)(x - x0)^2, (S_num - S_FL)/S_FL -> (13/72) l2/|l0| as
    q = |l0|/l2 -> infinity, with a residual of relative order 1/q.
 5. S_curv is NOT a bound: it lies above S_num for some potentials and below for others.
 6. Finiteness threshold of a variational estimate with trial shape g (phi = h g(r/R)) for the parabolic lambda:
    V < 0 for some h iff q > (1/2) Var(ln g) under the weight g^4 r^3 dr. For g_n = (1 + s^(2n))^(-1/n) this is
    13/72 at n = 1 and decreases monotonically to 3/32 as n -> infinity (g = min(1, s^-2)).
 7. The n = 1 member of the two-parameter family reproduces S_CGM; the n -> infinity member has the closed form
    S_inf = 4 pi^2 / (l2 (q - 3/32)) at its optimal ln h = x0 + 1/4.
 8. The exact bounce satisfies V < 0 (Derrick), hence (1/2) Var(ln phi) < q under phi^4 r^3 dr for the parabolic
    lambda: as q -> 0 the bounce must concentrate its phi^4 r^3 weight at one field value (thin wall).
Output: validation/results/derivation_checks.txt
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
from smvac.reference import analytic as A, bounce as B, potential as P, sm  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "results", "derivation_checks.txt")
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def radial_ratio(phi, dphi):
    """Q[phi] = int phi^4 d^4x / (int |grad phi|^2 d^4x)^2 for a radial profile (2 pi^2 r^3 dr measure)."""
    i4 = quad(lambda r: 2 * math.pi ** 2 * r ** 3 * phi(r) ** 4, 0, np.inf, limit=400)[0]
    i2 = quad(lambda r: 2 * math.pi ** 2 * r ** 3 * dphi(r) ** 2, 0, np.inf, limit=400)[0]
    return i4 / i2 ** 2


say("1. Sharp Sobolev constant K = 3/(32 pi^2) = %.15f" % (3 / (32 * math.pi ** 2)))
K = 3 / (32 * math.pi ** 2)
trials = {
    "FL 1/(1+r^2)": (lambda r: 1 / (1 + r * r), lambda r: -2 * r / (1 + r * r) ** 2),
    "(1+r^2)^-1.5": (lambda r: (1 + r * r) ** -1.5, lambda r: -3 * r * (1 + r * r) ** -2.5),
    "(1+r^2)^-3": (lambda r: (1 + r * r) ** -3, lambda r: -6 * r * (1 + r * r) ** -4),
    "Gaussian": (lambda r: math.exp(-r * r), lambda r: -2 * r * math.exp(-r * r)),
    "sech(r)": (lambda r: 2 * math.exp(-r) / (1 + math.exp(-2 * r)),                   # overflow-safe sech
                lambda r: -math.tanh(r) * 2 * math.exp(-r) / (1 + math.exp(-2 * r))),
    "(1+r) e^-r": (lambda r: (1 + r) * math.exp(-r), lambda r: -r * math.exp(-r)),
}
for name, (f, df) in trials.items():
    q = radial_ratio(f, df)
    say(f"   {name:14s} Q/K = {q / K:.12f}  {'<= 1 ok' if q <= K * (1 + 1e-10) else 'VIOLATION'}")

say()
say("2. Variational formula: direct r-space T^2/(4|V|) vs (8pi^2/3)/|lambdabar(h)|")
pots = {"parabola q=2.5": (P.parabola_model(-2.5 * 5e-5, 5e-5, 39.0), (25, 55)),
        "SM cw1 central": (P.sm_potential(sm.Inputs(125.15, 173.34), "cw1"), None)}
for name, (pot, search) in pots.items():
    an = A.estimates(pot, search)
    for dlnh in (-1.0, 0.0, 1.5):
        lnh = an.lnh_CGM + dlnh
        h = math.exp(lnh)
        T = quad(lambda r: 2 * math.pi ** 2 * r ** 3 * 0.5 * (2 * h * r / (1 + r * r) ** 2) ** 2, 0, np.inf, limit=400)[0]

        def vden(r):
            ph = h / (1 + r * r)
            return 2 * math.pi ** 2 * r ** 3 * 0.25 * pot.ev(math.log(ph))[0] * ph ** 4
        V = sum(quad(vden, a, b, limit=400, epsabs=0, epsrel=1e-13)[0]
                for a, b in [(0, 1), (1, 1e2), (1e2, 1e5), (1e5, 1e9), (1e9, 1e14)])
        direct = T * T / (4 * abs(V)) if V < 0 else float("inf")
        formula = A.EIGHT_PI2_3 / abs(A.lambda_bar(pot, lnh))
        say(f"   {name:15s} ln h = x_opt{dlnh:+.1f}: direct = {direct:.12g}  formula = {formula:.12g}  "
            f"rel diff = {abs(direct / formula - 1):.1e}")

say()
w_mean = quad(lambda w: 6 * w * (1 - w) * math.log(w), 0, 1)[0]
w_var = quad(lambda w: 6 * w * (1 - w) * math.log(w) ** 2, 0, 1)[0] - w_mean ** 2
say(f"3. Curvature coefficient: E[ln w] = {w_mean:.15f} (exact -5/6 = {-5/6:.15f}); "
    f"Var/2 = {w_var / 2:.15f} (exact 13/72 = {13/72:.15f})")

say()
say("4. First-order asymptotics for the parabolic model (numerical reference: M4, adaptive N)")
say("   q         (S-S_FL)/S_FL     13/(72 q)        ratio      (S_curv-S)/S   (S_CGM-S)/S")
for q in (3000.0, 1000.0, 300.0, 100.0, 30.0, 10.0, 3.0):
    l2 = 5e-5
    pot = P.parabola_model(-q * l2, l2, 39.0)
    an = A.estimates(pot, (20, 60))
    r = None
    for N in (256, 512, 1024):
        r = B.solve_m4(pot, an.lambda_min, an.lnh_CGM, N=N)
        if r.ok and abs(r.derrick) < 1e-10 and abs(r.dilatation) < 1e-10:
            break
    rel = (r.S - an.S_FL) / an.S_FL
    pred = 13 / (72 * q)
    say(f"   {q:7.0f}   {rel:.6e}    {pred:.6e}    {rel / pred:.6f}   {(an.S_curv - r.S) / r.S:+.3e}    "
        f"{(an.S_CGM - r.S) / r.S:+.3e}")

say()
say("5. S_curv is an approximation, not a bound: sign of (S_curv - S_num) at benchmark points")
for name, pot, search in [("parabola q=250", P.parabola_model(-0.0125, 5e-5, 39.0), (25, 55)),
                          ("SM cw1 central", P.sm_potential(sm.Inputs(125.15, 173.34), "cw1"), None)]:
    an = A.estimates(pot, search)
    r = B.solve_m4(pot, an.lambda_min, an.lnh_CGM, N=384)
    say(f"   {name:15s} S_FL={an.S_FL:.10g} S_num={r.S:.10g} S_curv={an.S_curv:.10g} S_CGM={an.S_CGM:.10g}  "
        f"S_curv {'>' if an.S_curv > r.S else '<'} S_num")

say()
say("6. Finiteness threshold q_c(n) = (1/2) Var(ln w_n) under w_n^4 s^3 ds, w_n = (1 + s^(2n))^(-1/n)")


def log_w(u, n):
    z = 2.0 * n * u
    return -(max(z, 0.0) + math.log1p(math.exp(-abs(z)))) / n


def qc(n):
    wt = lambda u: math.exp(4.0 * (u + log_w(u, n)))          # s^3 ds w^4 = s^4 w^4 du
    seg = ((-30, -1), (-1, 0), (0, 1), (1, 30))
    m = [sum(quad(lambda u, k=k: wt(u) * log_w(u, n) ** k, a, b, limit=500, epsabs=0, epsrel=1e-12)[0]
             for a, b in seg) for k in (0, 1, 2)]
    return 0.5 * (m[2] / m[0] - (m[1] / m[0]) ** 2)


prev = None
for n in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 1000.0, 1e4):
    v = qc(n)
    mono = "" if prev is None else ("  decreasing ok" if v < prev + 1e-12 else "  NOT MONOTONE")
    say(f"   n = {n:8g}   q_c = {v:.12f}{mono}")
    prev = v
say(f"   limits: 13/72 = {13 / 72:.12f}   3/32 = {3 / 32:.12f}")

say()
say("7. Two-parameter family: n = 1 member vs S_CGM; n -> infinity member vs 4 pi^2/(l2 (q - 3/32))")
l2 = 5e-5
for q in (2.5, 0.5, 0.15):
    pot = P.parabola_model(-q * l2, l2, 39.0)
    an = A.estimates(pot, (15, 65))
    if math.isfinite(an.S_CGM):
        T, V = A.cgm2_trial(pot, an.lnh_CGM, 1.0)
        say(f"   q = {q:5g}: n=1 member {T * T / (4 * abs(V)):.12g}   S_CGM {an.S_CGM:.12g}   rel diff "
            f"{abs(T * T / (4 * abs(V)) / an.S_CGM - 1):.1e}")
    T, V = A.cgm2_trial(pot, 39.25, 1e5)
    closed = 4 * math.pi ** 2 / (l2 * (q - 3 / 32))
    say(f"   q = {q:5g}: n=1e5 member at ln h = x0+1/4: {T * T / (4 * abs(V)):.10g}   closed form (n=inf) "
        f"{closed:.10g}   rel diff {abs(T * T / (4 * abs(V)) / closed - 1):.1e}")

say()
say("8. Concentration of the exact bounce: (1/2) Var(ln phi) under phi^4 r^3 dr vs q (parabolic model, M4)")
for q in (2.5, 0.5, 0.2, 0.1, 0.05):
    pot = P.parabola_model(-q * l2, l2, 39.0)
    an = A.estimates(pot, (15, 65))
    lnh = an.lnh_CGM if math.isfinite(an.lnh_CGM) else an.x_min + 5.0 / 6.0
    r = None
    for N in (256, 512, 1024):
        r = B.solve_m4(pot, an.lambda_min, lnh, N=N)
        if r.ok and abs(r.derrick) < 1e-9 and abs(r.dilatation) < 1e-9:
            break
    from scipy.interpolate import BarycentricInterpolator
    tn, fn = r.profile                                        # Chebyshev-Lobatto nodes: barycentric interpolation
    tau = np.linspace(tn.min(), tn.max(), 40001)
    f = BarycentricInterpolator(tn, fn)(tau)
    x = np.log(np.abs(f)) - tau                               # ln phi; the weight phi^4 r^3 dr = f^4 dtau
    wgt = np.full(tau.size, 2.0); wgt[1::2] = 4.0; wgt[0] = wgt[-1] = 1.0
    wgt *= (tau[1] - tau[0]) / 3.0 * f ** 4                   # Simpson weights times f^4
    mean = np.sum(wgt * x) / np.sum(wgt)
    half_var = 0.5 * (np.sum(wgt * (x - mean) ** 2) / np.sum(wgt))
    half_msd = 0.5 * (np.sum(wgt * (x - 39.0) ** 2) / np.sum(wgt))   # V < 0  <=>  half_msd < q (exact)
    say(f"   q = {q:5g}: (1/2) Var(ln phi) = {half_var:.6f}  (1/2) <(ln phi - x0)^2> = {half_msd:.6f}  "
        f"< q: {'ok' if half_msd < q else 'VIOLATION'}   ratio to q = {half_var / q:.4f}   "
        f"<ln phi> - x0 = {mean - 39.0:+.1e} (dilatation identity: 0)")

with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("wrote", OUT)
