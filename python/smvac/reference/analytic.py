"""Analytical approximations and bounds for the tree-level O(4) bounce action of V = lambda(phi) phi^4/4.

S_FL   = 8 pi^2 / (3 |lambda_min|)                         literature estimate; RIGOROUS LOWER BOUND
                                                           (Pohozaev identity + sharp Sobolev inequality)
S_curv = 8 pi^2 / (3 (|lambda_min| - 13/72 lambda''_min))  analytical approximation (first-order perturbation
                                                           theory about the FL instanton, quadratic lambda)
S_CGM  = 8 pi^2 / (3 |min_h lambdabar(h)|)                 variational upper estimate (a rigorous bound only if
                                                           the Coleman-Glaser-Martin infimum is attained)
with lambdabar(h) = 6 int_0^1 w (1 - w) lambda(h w) dw and x = ln(phi/GeV).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

EIGHT_PI2_3 = 8.0 * math.pi ** 2 / 3.0
CURV_COEFF = 13.0 / 72.0     # (1/2) Var(ln w) under the FL weight 6 w (1 - w) dw


def global_min(pot, search=None, n_grid: int = 4001):
    """Global minimum (lambda_min, x_min, at_boundary) of lambda(x) including frozen regions."""
    if getattr(pot, "x_lo", None) is not None:
        cands = [pot.x_lo, pot.x_hi]
        for r in pot.d1.roots():
            if abs(r.imag) < 1e-12 and pot.x_lo < r.real < pot.x_hi:
                cands.append(r.real)
        vals = [float(pot.lam(x)) for x in cands]
        i = int(np.argmin(vals))
        xm = cands[i]
        at_boundary = (abs(xm - pot.x_hi) < 1e-9) or (abs(xm - pot.x_lo) < 1e-9)
        return vals[i], xm, at_boundary
    a, b = search
    g = np.linspace(a, b, n_grid)
    v = pot.lam(g)
    i = int(np.argmin(v))
    r = minimize_scalar(lambda x: float(pot.lam(x)), bounds=(g[max(i - 1, 0)], g[min(i + 1, n_grid - 1)]),
                        method="bounded", options=dict(xatol=1e-12))
    return float(r.fun), float(r.x), i in (0, n_grid - 1)


def lambda_bar(pot, lnh: float) -> float:
    """lambdabar(h) = 6 int_{-inf}^0 e^{2y} (1 - e^y) lambda(ln h + y) dy (y = ln w), split at kinks."""
    f = lambda y: 6.0 * math.exp(2.0 * y) * (1.0 - math.exp(y)) * pot.ev(lnh + y)[0]
    pts = sorted(k - lnh for k in pot.kinks() if -45.0 < k - lnh < 0.0)
    edges = [-45.0] + pts + [0.0]
    return sum(quad(f, a, b, epsabs=1e-15, epsrel=1e-13, limit=200)[0] for a, b in zip(edges[:-1], edges[1:]))


@dataclass
class AnalyticResult:
    lambda_min: float
    x_min: float
    min_at_boundary: bool
    lambda2: float = float("nan")          # d^2 lambda / dx^2 at the minimum
    S_FL: float = float("inf")
    S_curv: float = float("nan")
    S_CGM: float = float("inf")
    lnh_CGM: float = float("nan")          # optimal ln phi(0) of the FL trial
    CGM_interior: bool = False
    lnLambdaB_lit: float = float("nan")    # Buttazzo identification: Lambda_B = argmin of lambda (mu = Lambda_B)
    lnLambdaB_FL: float = float("nan")     # 1/R of the optimal FL instanton: Lambda_B = h sqrt(|lambda|/8)
    notes: list = field(default_factory=list)


def estimates(pot, search=None) -> AnalyticResult:
    lmin, xm, atb = global_min(pot, search)
    res = AnalyticResult(lambda_min=lmin, x_min=xm, min_at_boundary=atb)
    if lmin >= 0.0:
        res.notes.append("absolutely stable: lambda_min >= 0")
        return res
    res.S_FL = EIGHT_PI2_3 / abs(lmin)
    res.lnLambdaB_lit = xm
    if not atb:
        res.lambda2 = float(pot.d2lam(xm))
        den = abs(lmin) - CURV_COEFF * res.lambda2
        res.S_curv = EIGHT_PI2_3 / den if den > 0 else float("inf")
        if den <= 0:
            res.notes.append("curvature correction exceeds |lambda_min|: first-order expansion invalid")
    else:
        res.notes.append("lambda minimum on a frozen boundary (Planck-sensitive / degenerate)")
    # CGM variational estimate: minimise lambdabar over ln h (scan, then bounded Brent)
    lo, hi = xm - 6.0, xm + 10.0
    grid = np.linspace(lo, hi, 161)
    vals = np.array([lambda_bar(pot, g) for g in grid])
    i = int(np.argmin(vals))
    if 0 < i < len(grid) - 1:
        r = minimize_scalar(lambda g: lambda_bar(pot, g), bounds=(grid[i - 1], grid[i + 1]), method="bounded",
                            options=dict(xatol=1e-9))
        lb, lnh = (r.fun, r.x) if r.fun <= vals[i] else (vals[i], grid[i])
        res.CGM_interior = True
    else:
        lb, lnh = vals[i], grid[i]
        res.notes.append("CGM optimum at the edge of the search range (no interior optimum; degenerate)")
    res.lnh_CGM = float(lnh)
    if lb < 0:
        res.S_CGM = EIGHT_PI2_3 / abs(lb)
        res.lnLambdaB_FL = float(lnh + 0.5 * math.log(abs(lb) / 8.0))
    return res


# --------------------------------------------------------------------------------------------------------------------
# Two-parameter variational estimate S_CGM2 (tested as a possible improvement in the thin-wall regime).
# Trial family phi(r) = h w(r/R), w(s) = (1 + s^(2n))^(-1/n). n = 1 is the Fubini-Lipatov shape; n -> infinity gives
# w = min(1, s^-2): a flat core joined continuously to the FL tail h (R/r)^2. The tail amplitude is therefore always
# tied to the core value h; the family cannot represent a thin wall across which the field drops below h R^2/r^2.
# Within the Coleman-Glaser-Martin characterisation every trial gives T^2/(4|V|) >= S_b, so the minimum over (h, n)
# is again a variational upper estimate (conditional on CGM attainment, like S_CGM), never larger than S_CGM.
# T(h, n) = (2 pi^2 h^2 / n) B(2 + 1/n, 1/n) exactly (R = 1; the ratio T^2/|V| is independent of R).
# Finiteness: for lambda = l2 [ (x - x0)^2 / 2 - q ], V < 0 for some h iff q > (1/2) Var(ln w) under the weight
# w^4 s^3 ds. This threshold is 13/72 for n = 1 and decreases monotonically to 3/32 as n -> infinity
# (docs/SCIENTIFIC_METHODOLOGY.md section 6.6a; validation/verify_derivations.py), so S_CGM2 is infinite for
# q <= 3/32: the second shape parameter moves the divergence of S_CGM but does not remove it.
# --------------------------------------------------------------------------------------------------------------------
def _log_w(u: float, n: float) -> float:
    """ln w(e^u) = -softplus(2 n u) / n, evaluated without overflow."""
    z = 2.0 * n * u
    return -(max(z, 0.0) + math.log1p(math.exp(-abs(z)))) / n


def _u_of_log_w(lw: float, n: float) -> float:
    """Inverse of _log_w for lw < 0."""
    L = -lw * n
    return (L + math.log(-math.expm1(-L))) / (2.0 * n)


def cgm2_trial(pot, lnh: float, n: float, epsrel: float = 1e-11):
    """(T, V) of the trial phi = h w(r) at R = 1 (units GeV^-1): S(R) = T R^2 + V R^4 is maximal at
    R^2 = T / (2|V|), where it equals T^2 / (4|V|) (V < 0 required)."""
    from scipy.special import beta as Bfun
    T = 2.0 * math.pi ** 2 * math.exp(2.0 * lnh) / n * Bfun(2.0 + 1.0 / n, 1.0 / n)

    def integrand(u):          # s^4 w^4 lambda(ln h + ln w), s = e^u  (V = (pi^2/2) h^4 int ... du)
        lw = _log_w(u, n)
        return math.exp(4.0 * (u + lw)) * pot.ev(lnh + lw)[0]
    pts = {-4.0 / n, 0.0, 4.0 / n}                                  # core-tail transition (width ~ 1/n)
    pts.update(_u_of_log_w(k - lnh, n) for k in pot.kinks() if k < lnh)   # images of kinks of lambda
    edges = [-18.0] + sorted(p for p in pts if -18.0 < p < 18.0) + [18.0]
    V = 0.5 * math.pi ** 2 * math.exp(4.0 * lnh) * sum(
        quad(integrand, a, b, limit=400, epsabs=0.0, epsrel=epsrel)[0] for a, b in zip(edges[:-1], edges[1:]))
    return T, V


def _cgm2_ratio(pot, lnh: float, n: float, epsrel: float = 1e-11) -> float:
    """T^2 / (4 |V|) for the trial (h = e^lnh, n); 1e300 if V >= 0 (finite penalty for Nelder-Mead)."""
    T, V = cgm2_trial(pot, lnh, n, epsrel)
    return T * T / (4.0 * abs(V)) if V < 0 else 1e300


def cgm2_estimate(pot, lnh_guess: float, n_max: float = 1000.0):
    """min over (ln h, n), 0.5 <= n <= n_max, of T^2/(4|V|) for the family above. Returns (S_CGM2, ln h, n).

    A coarse scan in ln h (lnh_guess - 2 ... + 1.5) and n (1 ... n_max) locates the region with V < 0; Nelder-Mead
    on ln S in (ln h, ln n) then refines the three best scan points. Returns (inf, lnh_guess, nan) if no member has
    V < 0. An optimum with n at n_max signals that the n -> infinity limit is approached (close to q = 3/32)."""
    from scipy.optimize import minimize
    if not math.isfinite(lnh_guess):
        return float("inf"), float("nan"), float("nan")
    ns = [n for n in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 300.0, 1000.0) if n <= n_max]
    grid = [(lnh_guess + d, n) for d in np.arange(-2.0, 1.5001, 0.1) for n in ns]
    vals = np.array([_cgm2_ratio(pot, l, n, epsrel=1e-8) for l, n in grid])
    if vals.min() >= 1e299:
        return float("inf"), float(lnh_guess), float("nan")
    clamp = lambda ln_n: min(n_max, max(0.5, math.exp(ln_n)))
    f = lambda p: math.log(_cgm2_ratio(pot, p[0], clamp(p[1])))
    best = (float("inf"), float("nan"), float("nan"))
    for k in np.argsort(vals)[:3]:
        if vals[k] >= 1e299:
            break
        l0, n0 = grid[k]
        r = minimize(f, [l0, math.log(n0)], method="Nelder-Mead",
                     options=dict(xatol=1e-8, fatol=1e-13, maxiter=4000,
                                  initial_simplex=[[l0, math.log(n0)], [l0 + 0.05, math.log(n0)],
                                                   [l0, math.log(n0) + 0.1]]))
        if math.exp(r.fun) < best[0]:
            best = (float(math.exp(r.fun)), float(r.x[0]), float(clamp(r.x[1])))
    return best
