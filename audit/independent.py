"""Independent reference calculations used to audit the legacy solver.

Independent of the legacy *numerics* (integrator, interpolation, quadrature,
optimiser); it re-uses only the legacy *physics inputs* (matching, beta functions,
lambda_eff) through ``audit.legacy_port`` so that numerical-method differences can be
isolated from physics differences.  Physics switches are passed through ``Options``.

Contents
--------
RGSolution         adaptive DOP853 solution of the RGEs with C^1-smooth dense output
                   (rtol ~1e-12), lambda_eff(t) and its exact t-derivative (chain rule).
tree_bounds        rigorous lower bound (Pohozaev + sharp Sobolev) and variational upper
                   estimate (rigorous only if the CGM infimum is attained) on the action of
                   V(phi) = lambda_eff(phi) phi^4 / 4:
                       8 pi^2/(3|min lambda|) <= S_b <= 8 pi^2/(3|min_h lambdabar(h)|),
                   lambdabar(h) = 6 int_0^1 w(1-w) lambda(h w) dw   (FL trial profile,
                   Coleman-Glaser-Martin variational characterisation).
shoot_bounce       numerical O(4) bounce by overshoot/undershoot shooting in the
                   scale-covariant variables f = r phi, tau = ln r.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar, brentq

from audit import legacy_port as L
from audit import _legacy_betas_extracted as B

PI2 = math.pi ** 2
EIGHT_PI2_3 = 8.0 * PI2 / 3.0


class RGSolution:
    """Dense, high-accuracy solution of the legacy RGEs from mu = M_t to mu = M_Pl."""

    def __init__(self, Mh: float, Mt: float, opt: L.Options = L.Options(), rtol: float = 1e-12,
                 atol: float = 1e-15, t_end: float | None = None, y0=None):
        self.opt = opt
        self.Mh, self.Mt = Mh, Mt
        self.t0 = 2.0 * math.log(Mt)
        self.tP = 2.0 * math.log(L.M_PLANCK) if t_end is None else t_end
        y0 = np.array(L.nnlo_matching(Mh, Mt, opt.matching) if y0 is None else y0, dtype=float)

        def rhs(t, y):
            return np.array(L.derivs(tuple(y), opt))

        def blowup(t, y):  # same perturbativity guard as legacy
            return min(4 * math.pi - abs(y[6]), 4 * math.pi - y[3])
        blowup.terminal = True
        self.sol = solve_ivp(rhs, (self.t0, self.tP), y0, method="DOP853", rtol=rtol, atol=atol,
                             dense_output=True, events=blowup)
        self.ok = self.sol.status == 0
        self.t_last = self.sol.t[-1]
        self.nfev = self.sol.nfev

    # couplings at t (vectorised), frozen outside [t0, t_last]
    def c(self, t):
        t = np.clip(t, self.t0, self.t_last)
        return self.sol.sol(t)

    def lambda_eff(self, t):
        cc = self.c(t)
        g1, g2, g3, yt, yb, ytau, lam = cc
        kz = 0.6 if self.opt.fix_z_coupling else 1.0
        xt, xw = yt * yt, g2 * g2
        xz = kz * g1 * g1 + xw
        corr = (-3.0 * xt**2 * (np.log(0.5 * xt) - 1.5) + 0.375 * xw**2 * (np.log(0.25 * xw) - 5.0 / 6.0)
                + 0.1875 * xz**2 * (np.log(0.25 * xz) - 5.0 / 6.0)) / (16.0 * PI2)
        return lam + corr

    def dlambda_eff_dt(self, t):
        """Exact d lambda_eff / d t (t = ln mu^2) by the chain rule; zero where frozen."""
        tt = np.atleast_1d(np.asarray(t, dtype=float))
        out = np.zeros_like(tt)
        inside = (tt > self.t0) & (tt < self.t_last)
        for i in np.nonzero(inside)[0]:
            cc = tuple(self.c(tt[i]))
            g1, g2, g3, yt, yb, ytau, lam = cc
            kz = 0.6 if self.opt.fix_z_coupling else 1.0
            xt, xw = yt * yt, g2 * g2
            xz = kz * g1 * g1 + xw
            dxt, dxw, dx1 = B.betaYt2(*cc), B.betaG2sq(*cc), B.betaG1sq(*cc)
            dxz = kz * dx1 + dxw
            dcorr = (-6.0 * xt * (math.log(0.5 * xt) - 1.0) * dxt
                     + 0.75 * xw * (math.log(0.25 * xw) - 1.0 / 3.0) * dxw
                     + 0.375 * xz * (math.log(0.25 * xz) - 1.0 / 3.0) * dxz) / (16.0 * PI2)
            out[i] = B.betaLambda(*cc) + dcorr
        return out if np.ndim(t) else float(out[0])


class FieldPotential:
    """lambda(phi) for V = lambda(phi) phi^4/4 with the legacy RG-improvement mu = phi.

    Below phi_low = max(M_t, v) and above M_Pl the couplings are frozen (legacy freezes
    above M_Pl; below v the legacy sets V = 0, which is irrelevant for the bounce)."""

    def __init__(self, rg: RGSolution, include_phi6: bool = False, phi6_scale: float = L.M_PLANCK,
                 phi6_coeff: float = 1.0):
        self.rg = rg
        self.t_lo = max(rg.t0, 2.0 * math.log(L.V_EW))
        self.t_hi = rg.t_last
        self.include_phi6 = include_phi6
        self.M = phi6_scale
        self.c6 = phi6_coeff

    def lam(self, lnphi):
        t = np.clip(2.0 * np.asarray(lnphi, dtype=float), self.t_lo, self.t_hi)
        val = self.rg.lambda_eff(t)
        if self.include_phi6:  # V6 = c6 phi^6/(6 M^2)  ->  lambda += (2/3) c6 phi^2/M^2
            val = val + (2.0 / 3.0) * self.c6 * np.exp(2.0 * np.asarray(lnphi) - 2.0 * math.log(self.M))
        return val

    def beta(self, lnphi):
        """d lambda / d ln phi."""
        lnphi = np.asarray(lnphi, dtype=float)
        t = 2.0 * lnphi
        inside = (t > self.t_lo) & (t < self.t_hi)
        b = np.where(inside, 2.0 * np.asarray(self.rg.dlambda_eff_dt(np.clip(t, self.t_lo, self.t_hi))), 0.0)
        if self.include_phi6:
            b = b + (4.0 / 3.0) * self.c6 * np.exp(2.0 * lnphi - 2.0 * math.log(self.M))
        return b


def lambda_min(pot: FieldPotential, n_grid: int = 4001):
    """Global minimum of lambda(phi) over the field range covered by the RG solution."""
    ln_lo, ln_hi = pot.t_lo / 2.0, max(pot.t_hi / 2.0, math.log(L.M_PLANCK)) + (6.0 if pot.include_phi6 else 0.0)
    grid = np.linspace(ln_lo, ln_hi, n_grid)
    vals = pot.lam(grid)
    i = int(np.argmin(vals))
    a, b = grid[max(i - 1, 0)], grid[min(i + 1, n_grid - 1)]
    if b > a:
        r = minimize_scalar(lambda x: float(pot.lam(x)), bounds=(a, b), method="bounded",
                            options=dict(xatol=1e-10))
        if r.fun <= vals[i]:
            return float(r.fun), float(r.x)
    return float(vals[i]), float(grid[i])


def lambda_bar(pot: FieldPotential, lnh: float, epsabs: float = 1e-13) -> float:
    """lambdabar(h) = 6 int_0^1 w(1-w) lambda(h w) dw, computed in x = ln w."""
    f = lambda x: 6.0 * math.exp(2 * x) * (1.0 - math.exp(x)) * float(pot.lam(lnh + x))
    val, _ = quad(f, -60.0, 0.0, epsabs=epsabs, epsrel=1e-12, limit=400)
    return val


def tree_bounds(pot: FieldPotential) -> dict:
    lmin, lnphi_min = lambda_min(pot)
    if lmin >= 0:
        return dict(lambda_min=lmin, S_lower=math.inf, S_upper=math.inf)
    # minimise lambdabar over ln h: scan then refine
    grid = np.linspace(lnphi_min - 6.0, lnphi_min + 8.0, 141)
    vals = np.array([lambda_bar(pot, x) for x in grid])
    i = int(np.argmin(vals))
    r = minimize_scalar(lambda x: lambda_bar(pot, x), bounds=(grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]),
                        method="bounded", options=dict(xatol=1e-8))
    lbar = min(r.fun, vals[i])
    return dict(lambda_min=lmin, lnphi_min=lnphi_min, lambdabar_min=lbar, lnh_opt=float(r.x),
                S_lower=EIGHT_PI2_3 / abs(lmin),
                S_upper=EIGHT_PI2_3 / abs(lbar) if lbar < 0 else math.inf)


# ----------------------------------------------------------------------------------------
# Numerical O(4) bounce by shooting.  Variables: tau = ln r, f = r phi, so phi = f e^{-tau}.
#   f'' = f + f^3 G(ln phi),  G = lambda + beta/4  (beta = d lambda / d ln phi)
# Mechanical energy E = f'^2/2 - f^2/2 - lambda f^4/4 obeys dE/dtau = beta f^4/4.
# Homoclinic orbit (bounce): E -> 0 at both ends.  Overshoot: f crosses 0.  Undershoot:
# f' returns to 0 at f > 0 after the turning point.
# ----------------------------------------------------------------------------------------
@dataclass
class ShotResult:
    kind: str           # 'over', 'under', 'undecided'
    sol: object
    tau_end: float


def _shoot(pot: FieldPotential, lnphi0: float, f_start: float, rtol: float, atol: float, tau_span: float = 80.0):
    phi0 = math.exp(lnphi0)
    G0 = float(pot.lam(lnphi0) + 0.25 * pot.beta(lnphi0))
    tau0 = math.log(f_start / phi0)
    r = math.exp(tau0)
    f_init = phi0 * r + G0 * phi0**3 * r**3 / 8.0
    fp_init = phi0 * r + 3.0 * G0 * phi0**3 * r**3 / 8.0

    def rhs(tau, y):
        f, fp = y
        lnphi = math.log(abs(f)) - tau if f != 0 else -1e300
        G = float(pot.lam(lnphi) + 0.25 * pot.beta(lnphi))
        return [fp, f + G * f**3]

    def ev_over(tau, y):
        return y[0]
    ev_over.terminal, ev_over.direction = True, -1

    def ev_under(tau, y):  # f' crosses zero upward (minimum of f) after the peak
        return y[1]
    ev_under.terminal, ev_under.direction = True, 1

    sol = solve_ivp(rhs, (tau0, tau0 + tau_span), [f_init, fp_init], method="DOP853", rtol=rtol, atol=atol,
                    dense_output=True, events=[ev_over, ev_under])
    if sol.t_events[0].size:
        return ShotResult("over", sol, sol.t_events[0][0])
    if sol.t_events[1].size:
        return ShotResult("under", sol, sol.t_events[1][0])
    return ShotResult("undecided", sol, sol.t[-1])


def _action_along(pot: FieldPotential, sol, tau_a: float, tau_b: float, n: int = 20001) -> dict:
    """S = T + V with T = 2pi^2 int (f'-f)^2/2, V = 2pi^2 int lambda f^4/4; plus virial checks."""
    taus = np.linspace(tau_a, tau_b, n)
    f, fp = sol.sol(taus)
    f = np.maximum(f, 1e-300)
    lnphi = np.log(f) - taus
    lam = pot.lam(lnphi)
    beta = pot.beta(lnphi)
    T = 2 * PI2 * np.trapezoid(0.5 * (fp - f) ** 2, taus)
    V = 2 * PI2 * np.trapezoid(0.25 * lam * f**4, taus)
    D = 2 * PI2 * np.trapezoid(0.25 * beta * f**4, taus)  # dilatation identity: should vanish
    return dict(S=T + V, T=T, V=V, derrick=(T + 2 * V) / T, dilatation=D / T)


def shoot_bounce(pot: FieldPotential, lnphi_lo: float, lnphi_hi: float, f_start: float = 1e-7,
                 rtol: float = 1e-11, atol: float = 1e-14, iters: int = 200) -> dict:
    """Bisection on ln phi(0) between an undershooting and an overshooting initial value."""
    lo, hi = _shoot(pot, lnphi_lo, f_start, rtol, atol), _shoot(pot, lnphi_hi, f_start, rtol, atol)
    if not (lo.kind == "under" and hi.kind == "over"):
        return dict(ok=False, reason=f"bracket not under/over: {lo.kind}/{hi.kind}")
    a, b = lnphi_lo, lnphi_hi
    best = None
    for _ in range(iters):
        m = 0.5 * (a + b)
        if m in (a, b):
            break
        s = _shoot(pot, m, f_start, rtol, atol)
        if s.kind == "under":
            a = m
        elif s.kind == "over":
            b = m
        else:
            break
        best = s
    # use the trajectory at the undershoot side of the boundary; integrate action up to the
    # point of closest approach to f = 0 (turning point / crossing)
    s = _shoot(pot, a, f_start, rtol, atol)
    sol = s.sol
    tau0 = sol.t[0]
    # cut where f has decayed to f_start-level relative to its max after the peak
    taus = np.linspace(tau0, s.tau_end, 200001)
    f = sol.sol(taus)[0]
    ipk = int(np.argmax(f))
    after = np.nonzero(f[ipk:] < f_start)[0]
    tau_cut = taus[ipk + after[0]] if after.size else s.tau_end
    act = _action_along(pot, sol, tau0, tau_cut)
    act.update(ok=True, lnphi0=a, phi0=math.exp(a), bracket_width=b - a, f_max=float(f[ipk]),
               lnphi_at_peak=float(math.log(f[ipk]) - taus[ipk]))
    return act


def find_shooting_bracket(pot: FieldPotential, ln_lo: float, ln_hi: float, n: int = 60, f_start: float = 1e-7,
                          rtol: float = 1e-10, atol: float = 1e-14):
    """Scan ln phi(0) and return all (under -> over) transitions."""
    grid = np.linspace(ln_lo, ln_hi, n)
    kinds = [_shoot(pot, x, f_start, rtol, atol).kind for x in grid]
    brackets = [(grid[i], grid[i + 1]) for i in range(n - 1) if kinds[i] == "under" and kinds[i + 1] == "over"]
    return brackets, list(zip(grid.tolist(), kinds))


class TabulatedPotential:
    """Fast lambda(ln phi) from a FieldPotential: uniform-grid cubic Hermite interpolant built
    from exact values and exact derivatives (beta).  lam/beta are evaluated from the SAME
    interpolant, so the bounce ODE force is exactly consistent with the potential and the
    Derrick/dilatation identities remain exact checks of the solver."""

    def __init__(self, pot: FieldPotential, ln_lo: float, ln_hi: float, h: float = 2e-3):
        self.include_phi6 = pot.include_phi6
        self.t_lo, self.t_hi = pot.t_lo, pot.t_hi
        n = int(math.ceil((ln_hi - ln_lo) / h))
        self.x0, self.h, self.n = ln_lo, (ln_hi - ln_lo) / n, n
        xs = ln_lo + self.h * np.arange(n + 1)
        self.y = [float(v) for v in pot.lam(xs)]
        self.d = [float(v) for v in pot.beta(xs)]

    def _eval(self, x):
        u = (x - self.x0) / self.h
        i = int(u)
        if i < 0:
            return self.y[0], 0.0
        if i >= self.n:
            return self.y[-1], 0.0
        s = u - i
        h = self.h
        y0, y1, d0, d1 = self.y[i], self.y[i + 1], self.d[i] * h, self.d[i + 1] * h
        s2, s3 = s * s, s * s * s
        val = (2 * s3 - 3 * s2 + 1) * y0 + (s3 - 2 * s2 + s) * d0 + (-2 * s3 + 3 * s2) * y1 + (s3 - s2) * d1
        der = ((6 * s2 - 6 * s) * y0 + (3 * s2 - 4 * s + 1) * d0 + (-6 * s2 + 6 * s) * y1 + (3 * s2 - 2 * s) * d1) / h
        return val, der

    def lam(self, x):
        if np.ndim(x):
            return np.array([self._eval(float(v))[0] for v in np.ravel(x)]).reshape(np.shape(x))
        return self._eval(float(x))[0]

    def beta(self, x):
        if np.ndim(x):
            return np.array([self._eval(float(v))[1] for v in np.ravel(x)]).reshape(np.shape(x))
        return self._eval(float(x))[1]
