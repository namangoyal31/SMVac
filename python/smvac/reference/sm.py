"""Standard Model inputs for the Python reference implementation.

Conventions (Buttazzo et al. 2013, arXiv:1307.3536, Sec. 3 and App. B):
  * V = -m^2 |H|^2 / 2 + lambda |H|^4, so M_h^2 = 2 lambda v^2 at tree level;
  * t = ln(mu^2); beta functions give d(x^2)/dt for gauge and Yukawa couplings and d(lambda)/dt;
  * g1 is GUT-normalised, g1^2 = (5/3) gY^2;
  * state vector c = (g1, g2, g3, yt, yb, ytau, lam).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from . import _betas_legacy_text as _B

PI2 = math.pi ** 2
LOOP1 = 16.0 * PI2
_B.LOOP1, _B.LOOP2, _B.LOOP3, _B.LOOP4 = LOOP1, LOOP1 ** 2, LOOP1 ** 3, LOOP1 ** 4

V_EW = 246.22             # GeV; v used in the tree-level relations (legacy value)
M_PLANCK = 1.22e19        # GeV; non-reduced Planck mass, upper end of the running
MB_TREE, MTAU_TREE = 4.0, 1.777   # GeV; tree-level y_b, y_tau inputs (legacy; effect < 1e-5 on lambda)


@dataclass(frozen=True)
class Inputs:
    Mh: float                 # Higgs pole mass [GeV]
    Mt: float                 # top pole mass [GeV]
    alpha3: float = 0.1184    # alpha_s(M_Z)
    MW: float = 80.384        # W pole mass [GeV]
    matching: str = "buttazzo2013"


def matching(inp: Inputs) -> tuple:
    """MS-bar couplings at mu = M_t."""
    if inp.matching in ("buttazzo2013", "extended"):   # eqs. (55), (57)-(60)
        das = (inp.alpha3 - 0.1184) / 0.0007
        dmw = (inp.MW - 80.384) / 0.014
        dt, dh = inp.Mt - 173.34, inp.Mh - 125.15
        lam = 0.12604 + 0.00206 * dh - 0.00004 * dt
        if inp.matching == "extended":   # tree-level M_h^2/(2 v^2) shape, same central value (src/matching.cpp)
            lam = inp.Mh ** 2 / (2 * V_EW ** 2) + (0.12604 - 125.15 ** 2 / (2 * V_EW ** 2)) - 0.00004 * dt
        yt = 0.93690 + 0.00556 * dt - 0.00042 * das
        g3 = 1.1666 + 0.00314 * das - 0.00046 * dt
        g2 = 0.64779 + 0.00004 * dt + 0.00011 * dmw
        gY = 0.35830 + 0.00011 * dt - 0.00020 * dmw
        g1 = math.sqrt(5.0 / 3.0) * gY
    elif inp.matching == "legacy":       # original code (docs/SCIENTIFIC_METHODOLOGY.md sec. 2)
        das = (inp.alpha3 - 0.1184) / 0.0007
        lam = inp.Mh ** 2 / (2 * V_EW ** 2) - 0.00313 - 0.00004 * (inp.Mt - 173.34)
        yt = math.sqrt(2.0) * inp.Mt / V_EW - 0.0587 - 0.00042 * das
        g1 = 0.46266 + 0.00006 * (inp.Mt - 173.34)
        g2 = 0.65355 + 0.00002 * (inp.Mt - 173.34)
        g3 = 1.1666 + 0.00314 * das - 0.00046 * (inp.Mt - 173.34)
    else:
        raise ValueError(f"unknown matching '{inp.matching}'")
    return (g1, g2, g3, yt, math.sqrt(2.0) * MB_TREE / V_EW, math.sqrt(2.0) * MTAU_TREE / V_EW, lam)


def beta_squares(c) -> tuple:
    """(d g1^2/dt, d g2^2/dt, d g3^2/dt, d yt^2/dt, d yb^2/dt, d ytau^2/dt, d lam/dt) — Buttazzo App. B."""
    g1, g2, g3, yt, yb, ytau, lam = c
    yt2 = yt * yt
    b_g3 = _B.betaG3sq(*c) - 2.0 * g3 ** 10 / _B.LOOP4 * 2472.28                  # 4-loop sign (eq. 98)
    b_lam = _B.betaLambda(*c) + (yt2 ** 3 - yt2 ** 2) * (250.494 * g3 * g3 + 74.138 * g2 * g2
                                                         + 33.930 * g1 * g1) / _B.LOOP3     # eq. 99
    b_yt = _B.betaYt2(*c) + yt2 / _B.LOOP3 * (yt2 ** 2 - yt2) * (
        58.6028 * yt2 + 198 * lam - 157 * g3 * g3 - 1593 * g2 * g2 / 16 - 2437 * g1 * g1 / 80)  # eq. 101
    return (_B.betaG1sq(*c), _B.betaG2sq(*c), b_g3, b_yt, _B.betaYb2(*c), _B.betaYtau2(*c), b_lam)


def derivs(c) -> list:
    b = beta_squares(c)
    return [b[i] / (2.0 * c[i]) for i in range(6)] + [b[6]]


def lambda_cw1(c):
    """One-loop Coleman-Weinberg effective quartic at mu = h (Landau gauge, MS-bar, h >> v),
    Buttazzo et al. eq. (106) with Gamma = 0 and without the Goldstone/Higgs term."""
    import numpy as np
    g1, g2, g3, yt, yb, ytau, lam = c
    xt, xw = yt * yt, g2 * g2
    xz = 0.6 * g1 * g1 + xw                       # g2^2 + gY^2
    return lam + (-3.0 * xt ** 2 * (np.log(0.5 * xt) - 1.5) + 0.375 * xw ** 2 * (np.log(0.25 * xw) - 5.0 / 6.0)
                  + 0.1875 * xz ** 2 * (np.log(0.25 * xz) - 5.0 / 6.0)) / LOOP1
