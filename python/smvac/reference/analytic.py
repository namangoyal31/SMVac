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
