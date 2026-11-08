"""Faithful Python port of the legacy solver in ``apps/solver_numerical.cpp``.

Purpose: audit tooling only.  It reproduces the legacy C++ computation exactly
(validated against the compiled legacy binary by ``audit/validate_port.py``) and
exposes switches that isolate individual suspected defects, so that the
numerical impact of each defect can be measured on its own.

The beta-function bodies are extracted mechanically from the C++ source into
``_legacy_betas_extracted.py``; they are NOT re-typed by hand.

Every switch defaults to the legacy behaviour.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from bisect import bisect_right

from . import _legacy_betas_extracted as _B

PI = math.pi
PI2 = PI * PI
LOOP1 = 16.0 * PI2
_B.LOOP1, _B.LOOP2, _B.LOOP3, _B.LOOP4 = LOOP1, LOOP1**2, LOOP1**3, LOOP1**4

V_EW = 246.22          # GeV, legacy value
M_PLANCK = 1.22e19     # GeV, legacy value (non-reduced Planck mass)
MTAU, MB = 1.777, 4.0  # GeV, legacy values
ALPHA3_MZ = 0.1184


@dataclass
class Options:
    fix_z_coupling: bool = False     # use (3/5) g1^2 + g2^2 for the Z mass (GUT-normalised g1)
    flip_4loop_sign: bool = False    # correct sign of the 4-loop QCD term
    include_phi6: bool = True        # legacy adds V6 = h^6 / (6 M_Pl^2) in the bounce integrand
    kinetic_cutoff: float = 5e5      # legacy early exit returns T (not S) above this
    n_quad: int = 2048
    gss_tol: float = 1e-13
    gss_half_width: float = 8.0      # bracket half-width in ln R
    dt: float = 0.1                  # RK4 step in t = ln mu^2
    interp_all: bool = False         # legacy interpolates only lambda; True interpolates all couplings
    fix_3loop_lambda: bool = False   # restore missing y_t^2 factor in 3 three-loop beta_lambda terms (BDGGSSS 2013 eq. 99)
    fix_3loop_yt: bool = False       # restore missing y_t^2 factor in the 3-loop y_t^4(...) group (eq. 101)
    matching: str = "legacy"         # "legacy" or "buttazzo2013" (eqs. 55, 57-60 of arXiv:1307.3536)


Couplings = tuple  # (g1, g2, g3, yt, yb, ytau, lam)


def nnlo_matching(Mh: float, Mt: float, matching: str = "legacy") -> Couplings:
    if matching == "buttazzo2013":
        d_as = (ALPHA3_MZ - 0.1184) / 0.0007
        dt, dh = Mt - 173.34, Mh - 125.15
        gY = 0.35830 + 0.00011 * dt
        return (math.sqrt(5.0 / 3.0) * gY, 0.64779 + 0.00004 * dt, 1.1666 + 0.00314 * d_as - 0.00046 * dt,
                0.93690 + 0.00556 * dt - 0.00042 * d_as, math.sqrt(2.0) * MB / V_EW, math.sqrt(2.0) * MTAU / V_EW,
                0.12604 + 0.00206 * dh - 0.00004 * dt)
    dAlphas = (ALPHA3_MZ - 0.1184) / 0.0007
    lambda_tree = (Mh * Mh) / (2.0 * V_EW * V_EW)
    yt_tree = (math.sqrt(2.0) * Mt) / V_EW
    delta_lambda = -0.00313 - 0.00004 * (Mt - 173.34)
    delta_yt = -0.0587 - 0.00042 * dAlphas
    g1 = 0.46266 + 0.00006 * (Mt - 173.34)
    g2 = 0.65355 + 0.00002 * (Mt - 173.34)
    g3 = 1.1666 + 0.00314 * dAlphas - 0.00046 * (Mt - 173.34)
    return (g1, g2, g3, yt_tree + delta_yt, math.sqrt(2.0) * MB / V_EW,
            math.sqrt(2.0) * MTAU / V_EW, lambda_tree + delta_lambda)


def lambda_eff(c: Couplings, opt: Options = Options()) -> float:
    g1, g2, g3, yt, yb, ytau, lam = c
    yt2 = yt * yt
    yt4 = yt2 * yt2
    term_t = -3.0 * yt4 * (math.log(0.5 * yt2) - 1.5)
    g2_2 = g2 * g2
    g2_4 = g2_2 * g2_2
    term_W = 0.375 * g2_4 * (math.log(0.25 * g2_2) - 5.0 / 6.0)
    g12 = (0.6 * g1 * g1 if opt.fix_z_coupling else g1 * g1) + g2_2
    g12_2 = g12 * g12
    term_Z = 0.1875 * g12_2 * (math.log(0.25 * g12) - 5.0 / 6.0)
    return lam + (term_t + term_W + term_Z) / (16.0 * PI2)


def derivs(c: Couplings, opt: Options = Options()) -> Couplings:
    g1, g2, g3, yt, yb, ytau, lam = c
    b3 = _B.betaG3sq(*c)
    if opt.flip_4loop_sign:
        b3 -= 2.0 * (g3**10) / _B.LOOP4 * 2472.28
    bl = _B.betaLambda(*c)
    byt = _B.betaYt2(*c)
    yt2 = yt * yt
    if opt.fix_3loop_lambda:
        bl += (yt2**3 - yt2**2) * (250.494 * g3 * g3 + 74.138 * g2 * g2 + 33.930 * g1 * g1) / _B.LOOP3
    if opt.fix_3loop_yt:
        byt += yt2 / _B.LOOP3 * (yt2**2 - yt2) * (58.6028 * yt2 + 198 * lam - 157 * g3 * g3
                                                   - 1593 * g2 * g2 / 16 - 2437 * g1 * g1 / 80)
    return (0.5 / g1 * _B.betaG1sq(*c), 0.5 / g2 * _B.betaG2sq(*c), 0.5 / g3 * b3,
            0.5 / yt * byt, 0.5 / yb * _B.betaYb2(*c), 0.5 / ytau * _B.betaYtau2(*c), bl)


def rk4_step(y: Couplings, dt: float, opt: Options = Options()) -> Couplings:
    k1 = derivs(y, opt)
    k2 = derivs(tuple(a + 0.5 * dt * b for a, b in zip(y, k1)), opt)
    k3 = derivs(tuple(a + 0.5 * dt * b for a, b in zip(y, k2)), opt)
    k4 = derivs(tuple(a + dt * b for a, b in zip(y, k3)), opt)
    return tuple(a + dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(y, k1, k2, k3, k4))


@dataclass
class RGETable:
    t: list = field(default_factory=list)
    y: list = field(default_factory=list)
    interp_all: bool = False

    def add(self, t: float, y: Couplings) -> None:
        self.t.append(t)
        self.y.append(y)

    def get(self, t: float) -> Couplings:
        """Legacy RGEHelper::get_params: linear interpolation of lambda only,
        all other couplings taken from the left grid point (piecewise constant)."""
        i = bisect_right(self.t, t)
        if i == 0:
            return self.y[0]
        if i == len(self.t):
            return self.y[-1]
        t1, t2 = self.t[i - 1], self.t[i]
        p1, p2 = self.y[i - 1], self.y[i]
        f = (t - t1) / (t2 - t1)
        if self.interp_all:
            return tuple(a + f * (b - a) for a, b in zip(p1, p2))
        return p1[:6] + (p1[6] + f * (p2[6] - p1[6]),)


@dataclass
class RunResult:
    status: int                 # 1 stable, 2 metastable, 3 unstable, 4 non-perturbative
    table: RGETable
    mu_inst: float
    t_min_lambda: float
    min_lambda_eff: float


def run_rge(Mh: float, Mt: float, opt: Options = Options()) -> RunResult:
    y = nnlo_matching(Mh, Mt, opt.matching)
    t0 = 2 * math.log(Mt)
    tP = 2 * math.log(M_PLANCK)
    tab = RGETable(interp_all=opt.interp_all)
    is_unstable = False
    mu_inst = M_PLANCK
    min_le = 1.0
    t_min = t0
    t = t0
    dt = opt.dt
    tab.add(t, y)
    t_v = 2 * math.log(V_EW)
    if t >= t_v:
        min_le = lambda_eff(y, opt)
        if min_le <= 0.0:
            mu_inst = math.exp(t / 2.0)
            is_unstable = True
    while t < tP:
        if t + dt > tP:
            dt = tP - t
        y = rk4_step(y, dt, opt)
        t += dt
        tab.add(t, y)
        if abs(y[6]) > 4 * PI or y[3] > 4 * PI or not (math.isfinite(y[6]) and math.isfinite(y[3])):
            return RunResult(4, tab, mu_inst, t_min, min_le)
        if t >= t_v:
            le = lambda_eff(y, opt)
            if le < min_le:
                min_le, t_min = le, t
            if not is_unstable and le < 0:
                mu_inst = math.exp(t / 2.0)
                is_unstable = True
    return RunResult(2 if is_unstable else 1, tab, mu_inst, t_min, min_le)


def _lambda_at_field(tab: RGETable, h: float, opt: Options) -> float:
    h_rge = min(h, M_PLANCK)
    if h_rge <= 246.22:
        return 0.0
    return lambda_eff(tab.get(2.0 * math.log(h_rge)), opt)


def action_at_R(tab: RGETable, mu_inst: float, R: float, opt: Options = Options(),
                return_parts: bool = False):
    """Legacy evaluate_action_at_R (the quantity the legacy code calls S_exact before minimisation)."""
    if R <= 0.0:
        return 1e100
    mu_R = mu_inst / R
    lam_R = tab.get(2.0 * math.log(mu_R))[6]
    if lam_R >= 0.0:
        return 1e100
    abs_lam = abs(lam_R)
    pref = math.sqrt(2.0 / abs_lam)
    kinetic = 16.0 * PI2 / (3.0 * abs_lam)
    if kinetic > opt.kinetic_cutoff:
        return kinetic
    N = opt.n_quad
    du = 2.0 / N
    pot = 0.0
    alpha = 4.0
    for i in range(1, N):
        u = -1.0 + i * du
        ratio = (1.0 + u) / (1.0 - u)
        e2x = ratio**alpha
        e4x = e2x * e2x
        phi = pref * 2.0 / (R * (e2x + 1.0))
        lam_phi = _lambda_at_field(tab, phi * mu_inst, opt)
        phi4 = phi**4
        V = 0.25 * lam_phi * phi4
        if opt.include_phi6:
            V += (1.0 / 6.0) * (phi * mu_inst / M_PLANCK) ** 2 * phi4
        f = 2.0 * PI2 * R**4 * e4x * V * alpha / (1.0 - u * u)
        pot += (4.0 / 3.0 if i % 2 == 1 else 2.0 / 3.0) * f * du
    if return_parts:
        return kinetic + pot, kinetic, pot
    return kinetic + pot


def golden_section(f, a: float, b: float, tol: float) -> tuple[float, float]:
    invphi = (math.sqrt(5.0) - 1.0) / 2.0
    invphi2 = (3.0 - math.sqrt(5.0)) / 2.0
    h = b - a
    if h <= tol:
        x = 0.5 * (a + b)
        return x, f(x)
    n = int(math.ceil(math.log(tol / h) / math.log(invphi)))
    c, d = a + invphi2 * h, a + invphi * h
    yc, yd = f(c), f(d)
    for _ in range(n):
        if yc < yd:
            b, d, yd = d, c, yc
            h = invphi * h
            c = a + invphi2 * h
            yc = f(c)
        else:
            a, c, yc = c, d, yd
            h = invphi * h
            d = a + invphi * h
            yd = f(d)
    x = 0.5 * (a + b)
    return x, f(x)


def threshold_legacy(mu_inst: float) -> float:
    Tv = 1.179e44 / V_EW
    return 4.0 * math.log(Tv * mu_inst)


def classify(Mh: float, Mt: float, opt: Options = Options()) -> dict:
    run = run_rge(Mh, Mt, opt)
    out = dict(Mh=Mh, Mt=Mt, status=run.status, S_exact=-1.0, S_approx=-1.0,
               mu_inst=run.mu_inst, mu_min=math.exp(run.t_min_lambda / 2), lambda_eff_min=run.min_lambda_eff)
    if run.status in (1, 4):
        return out
    S_approx = 8.0 * PI2 / (3.0 * abs(run.min_lambda_eff))
    mu_opt = math.exp(run.t_min_lambda / 2.0)
    logR = math.log(run.mu_inst / mu_opt)
    xmin, S_exact = golden_section(lambda lr: action_at_R(run.table, run.mu_inst, math.exp(lr), opt),
                                   logR - opt.gss_half_width, logR + opt.gss_half_width, opt.gss_tol)
    thr = threshold_legacy(run.mu_inst)
    out.update(S_exact=S_exact, S_approx=S_approx, S_threshold=thr, logR_min=xmin,
               mu_R_min=run.mu_inst / math.exp(xmin), status=2 if S_exact > thr else 3, run=run)
    return out
