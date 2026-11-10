"""RG-improved quartic lambda(x), x = ln(phi/GeV), and model potentials.

Every potential exposes the same interface, used by the analytic estimates and all bounce solvers:
    lam(x), dlam(x), d2lam(x)   vectorised: lambda, d lambda/dx, d^2 lambda/dx^2
    ev(x) -> (lambda, dlambda/dx)  fast scalar evaluation (pure Python)
    x_lo, x_hi                  couplings are frozen outside [x_lo, x_hi] (None for models)
"""
from __future__ import annotations

import math

import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.integrate import solve_ivp

from . import sm


class HermiteTable5:
    """Uniform-grid quintic Hermite interpolant (C^2) built from values and two derivatives."""

    def __init__(self, x0: float, x1: float, h: float, f, df, d2f):
        n = max(2, int(math.ceil((x1 - x0) / h)))
        self.x0, self.h, self.n = x0, (x1 - x0) / n, n
        xs = x0 + self.h * np.arange(n + 1)
        self.y = [float(v) for v in f(xs)]
        self.d = [float(v) * self.h for v in df(xs)]
        self.dd = [float(v) * self.h * self.h for v in d2f(xs)]

    def ev(self, x: float):
        u = (x - self.x0) / self.h
        i = int(u)
        if i < 0:
            i, s = 0, 0.0
        elif i >= self.n:
            i, s = self.n - 1, 1.0
        else:
            s = u - i
        s2 = s * s; s3 = s2 * s; s4 = s3 * s; s5 = s4 * s
        y0, y1, d0, d1, e0, e1 = self.y[i], self.y[i + 1], self.d[i], self.d[i + 1], self.dd[i], self.dd[i + 1]
        val = (y0 * (1 - 10 * s3 + 15 * s4 - 6 * s5) + d0 * (s - 6 * s3 + 8 * s4 - 3 * s5)
               + e0 * 0.5 * (s2 - 3 * s3 + 3 * s4 - s5) + y1 * (10 * s3 - 15 * s4 + 6 * s5)
               + d1 * (-4 * s3 + 7 * s4 - 3 * s5) + e1 * 0.5 * (s3 - 2 * s4 + s5))
        der = (y0 * (-30 * s2 + 60 * s3 - 30 * s4) + d0 * (1 - 18 * s2 + 32 * s3 - 15 * s4)
               + e0 * 0.5 * (2 * s - 9 * s2 + 12 * s3 - 5 * s4) + y1 * (30 * s2 - 60 * s3 + 30 * s4)
               + d1 * (-12 * s2 + 28 * s3 - 15 * s4) + e1 * 0.5 * (3 * s2 - 8 * s3 + 5 * s4)) / self.h
        return val, der


class FrozenChebyshevPotential:
    """lambda(x) as a Chebyshev series on [x_lo, x_hi]; constant (frozen couplings) outside."""

    def __init__(self, func, x_lo: float, x_hi: float, deg: int | None = None, table_h: float = 5e-3):
        self.x_lo, self.x_hi = x_lo, x_hi
        dom = [x_lo, x_hi]
        if deg is None:  # increase the degree until the tail coefficients are at round-off level
            for deg in (48, 64, 96, 128, 192, 256, 384):
                ser = C.Chebyshev.interpolate(func, deg, domain=dom)
                scale = np.max(np.abs(ser.coef))
                if np.max(np.abs(ser.coef[-6:])) < 1e-15 * scale:
                    break
        else:
            ser = C.Chebyshev.interpolate(func, deg, domain=dom)
        self.ser, self.deg = ser, deg
        self.tail = float(np.max(np.abs(ser.coef[-6:])) / np.max(np.abs(ser.coef)))
        self.d1, self.d2 = ser.deriv(1), ser.deriv(2)
        self.lam_lo, self.lam_hi = float(ser(x_lo)), float(ser(x_hi))
        pad = 80.0
        self._tab = HermiteTable5(x_lo, x_hi, table_h, self.ser, self.d1, self.d2)
        self._pad = pad

    def lam(self, x):
        x = np.asarray(x, dtype=float)
        return np.where(x < self.x_lo, self.lam_lo, np.where(x > self.x_hi, self.lam_hi,
                                                              self.ser(np.clip(x, self.x_lo, self.x_hi))))

    def dlam(self, x):
        x = np.asarray(x, dtype=float)
        inside = (x > self.x_lo) & (x < self.x_hi)
        return np.where(inside, self.d1(np.clip(x, self.x_lo, self.x_hi)), 0.0)

    def d2lam(self, x):
        x = np.asarray(x, dtype=float)
        inside = (x > self.x_lo) & (x < self.x_hi)
        return np.where(inside, self.d2(np.clip(x, self.x_lo, self.x_hi)), 0.0)

    def ev(self, x: float):
        if x <= self.x_lo:
            return self.lam_lo, 0.0
        if x >= self.x_hi:
            return self.lam_hi, 0.0
        return self._tab.ev(x)

    def kinks(self):
        return [self.x_lo, self.x_hi]


class ModelPotential:
    """Analytic model lambda(x) (no freezing) given callables for lambda and its derivatives."""

    def __init__(self, lam, dlam, d2lam, name: str = "model"):
        self._l, self._d, self._d2, self.name = lam, dlam, d2lam, name
        self.x_lo = self.x_hi = None

    def lam(self, x):
        return self._l(np.asarray(x, dtype=float))

    def dlam(self, x):
        return self._d(np.asarray(x, dtype=float))

    def d2lam(self, x):
        return self._d2(np.asarray(x, dtype=float))

    def ev(self, x: float):
        return float(self._l(x)), float(self._d(x))

    def kinks(self):
        return []


def constant_model(l0: float) -> ModelPotential:
    return ModelPotential(lambda x: l0 + 0 * x, lambda x: 0 * x, lambda x: 0 * x, f"const({l0})")


def parabola_model(l0: float, l2: float, x0: float) -> ModelPotential:
    """lambda = l0 + (l2/2)(x - x0)^2 : minimum l0 at x0 with curvature l2 (thick wall for |l0|/l2 >> 1,
    thin wall with a nearly degenerate second vacuum for |l0|/l2 -> 0)."""
    return ModelPotential(lambda x: l0 + 0.5 * l2 * (x - x0) ** 2, lambda x: l2 * (x - x0), lambda x: l2 + 0 * x,
                          f"parabola({l0},{l2},{x0})")


class RGSolution:
    """Dense DOP853 solution of the SM RGEs from mu = M_t up to M_Pl (or a perturbativity failure)."""

    def __init__(self, inp: sm.Inputs, rtol: float = 1e-13, atol: float = 1e-16):
        self.inp = inp
        self.t0 = 2.0 * math.log(inp.Mt)
        self.tP = 2.0 * math.log(sm.M_PLANCK)
        y0 = np.array(sm.matching(inp))

        def blowup(t, y):
            return min(4 * math.pi - abs(y[6]), 4 * math.pi - y[3])
        blowup.terminal = True
        self.sol = solve_ivp(lambda t, y: sm.derivs(y), (self.t0, self.tP), y0, method="DOP853",
                             rtol=rtol, atol=atol, dense_output=True, events=blowup)
        self.perturbative = self.sol.status == 0
        self.t_end = float(self.sol.t[-1])

    def c(self, t):
        return self.sol.sol(np.clip(t, self.t0, self.t_end))


def sm_potential(inp: sm.Inputs, prescription: str = "cw1", rg: RGSolution | None = None, **kw):
    """lambda_P(x) for prescription 'tree' (MS-bar lambda at mu = phi) or 'cw1' (one-loop CW at mu = phi)."""
    rg = rg or RGSolution(inp)
    x_lo = 0.5 * max(rg.t0, 2.0 * math.log(sm.V_EW))
    x_hi = 0.5 * rg.t_end
    if prescription == "tree":
        func = lambda x: rg.c(2.0 * x)[6]
    elif prescription == "cw1":
        func = lambda x: sm.lambda_cw1(rg.c(2.0 * x))
    else:
        raise ValueError(prescription)
    pot = FrozenChebyshevPotential(func, x_lo, x_hi, **kw)
    pot.rg, pot.prescription, pot.inp = rg, prescription, inp
    return pot
