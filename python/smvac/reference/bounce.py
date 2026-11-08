"""Candidate numerical methods for the O(4) bounce of V = lambda(phi) phi^4 / 4 (problem P).

Variables: tau = ln r, f = r phi, x = ln phi = ln f - tau. Bounce equation
    f'' = f + f^3 G(x),   G = lambda + beta/4,   beta = d lambda/dx,
with f -> 0 like e^{+tau} (r -> 0) and e^{-tau} (r -> infinity). Action
    S = 2 pi^2 int dtau [ (f' - f)^2/2 + lambda f^4/4 ] = T + V.
Checks for a solution: Derrick identity T + 2V = 0; dilatation identity int beta f^4 dtau = 0.

M0  legacy construction (first-order FL estimate with lambda_c = lambda(x = -ln R)), for comparison only
M1  overshoot/undershoot shooting in physical variables (s = r phi0, u = phi/phi0)
M2  overshoot/undershoot shooting in scale-covariant variables (tau, f)
M3  global collocation BVP (scipy.solve_bvp, 4th-order Lobatto IIIA, adaptive mesh), Robin end conditions
M4  Chebyshev spectral collocation with damped Newton, Robin end conditions
None of these is "exact": each returns a numerical solution whose accuracy is assessed separately.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
from scipy.optimize import minimize_scalar

TWO_PI2 = 2.0 * math.pi ** 2


@dataclass
class BounceResult:
    method: str
    ok: bool
    S: float = float("nan")
    T: float = float("nan")
    V: float = float("nan")
    derrick: float = float("nan")      # (T + 2V)/T
    dilatation: float = float("nan")   # int beta f^4 / int |lambda| f^4
    lnphi0: float = float("nan")
    lnR: float = float("nan")          # tau at the maximum of f = r phi (R for a Fubini-Lipatov profile)
    seconds: float = float("nan")
    info: dict = field(default_factory=dict)
    profile: tuple | None = None       # (tau, f) on the solver's nodes, when available


def _simpson(y, h):
    n = len(y) - 1
    if n % 2:
        raise ValueError("Simpson needs an even number of intervals")
    return h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum())


def _functionals(pot, tau, f, fp, weights):
    """T, V, S, identity residuals from a sampled profile (weights: quadrature weights on tau)."""
    fpos = np.maximum(f, 1e-300)
    x = np.log(fpos) - tau
    lam, beta = pot.lam(x), pot.dlam(x)
    T = TWO_PI2 * np.dot(weights, 0.5 * (fp - f) ** 2) + TWO_PI2 * f[-1] ** 2   # + right-tail f ~ e^{-tau}
    V = TWO_PI2 * np.dot(weights, 0.25 * lam * f ** 4)
    dil = np.dot(weights, beta * f ** 4) / np.dot(weights, np.abs(lam) * f ** 4)
    imax = int(np.argmax(f))
    return dict(T=T, V=V, S=T + V, derrick=(T + 2 * V) / T, dilatation=dil, lnR=float(tau[imax]),
                f_max=float(f[imax]))


def admissible_interval(pot, x_center, half=10.0, step=0.01):
    """Connected interval around x_center (the lambda minimum) where G = lambda + beta/4 < 0, i.e. V'(phi) < 0:
    the only starting values phi(0) from which a bounce can roll towards the false vacuum. Its upper end is the true
    vacuum (or the edge of the scan for potentials unbounded below)."""
    import numpy as _np
    xs = _np.arange(x_center - half, x_center + half + step, step)
    G = pot.lam(xs) + 0.25 * pot.dlam(xs)
    i0 = int(_np.argmin(_np.abs(xs - x_center)))
    if G[i0] >= 0:
        neg = _np.nonzero(G < 0)[0]
        if not neg.size:
            return None
        i0 = int(neg[_np.argmin(_np.abs(neg - i0))])
    lo = i0
    while lo > 0 and G[lo - 1] < 0:
        lo -= 1
    hi = i0
    while hi < len(xs) - 1 and G[hi + 1] < 0:
        hi += 1
    return float(xs[lo]), float(xs[hi])


# ----------------------------------------------------------------------------- M2: (tau, f) shooting
def _m2_shot(pot, x0, eps_f, rtol, atol, span, f_big=1e8):
    phi0 = math.exp(x0)
    l0, d0 = pot.ev(x0)
    G0 = l0 + 0.25 * d0
    tau0 = math.log(eps_f / phi0)
    r = eps_f / phi0
    y0 = [phi0 * r + G0 * phi0 ** 3 * r ** 3 / 8.0, phi0 * r + 3.0 * G0 * phi0 ** 3 * r ** 3 / 8.0]
    ev = pot.ev

    def rhs(t, y):
        f, fp = y
        if f <= 0.0:
            return [fp, f]
        l, d = ev(math.log(f) - t)
        return [fp, f + f * f * f * (l + 0.25 * d)]

    def over(t, y):
        return y[0]
    over.terminal, over.direction = True, -1

    def under(t, y):
        return y[1]
    under.terminal, under.direction = True, 1

    def escape(t, y):  # phi grows beyond phi0: started on the wrong side of the true vacuum
        return (math.log(max(y[0], 1e-300)) - t) - (x0 + 0.5)
    escape.terminal, escape.direction = True, 1

    def runaway(t, y):  # f = r phi grows exponentially while phi lingers near a true vacuum: thin-wall overshoot
        return y[0] - f_big
    runaway.terminal, runaway.direction = True, 1

    def turn(t, y):  # phi' = (f' - f)/r^2 becomes positive: the field turns back (undershoot)
        return (y[1] - y[0]) - 1e-9 * abs(y[0])   # threshold above round-off of f' - f near the start
    turn.terminal, turn.direction = True, 1

    sol = solve_ivp(rhs, (tau0, tau0 + span), y0, method="DOP853", rtol=rtol, atol=atol, dense_output=True,
                    events=[over, under, escape, runaway, turn])
    if sol.t_events[2].size:
        return "escape", sol, sol.t_events[2][0]
    if sol.t_events[3].size:
        return "over", sol, sol.t_events[3][0]
    if sol.t_events[4].size:
        return "under", sol, sol.t_events[4][0]
    if sol.t_events[0].size:
        return "over", sol, sol.t_events[0][0]
    if sol.t_events[1].size:
        return "under", sol, sol.t_events[1][0]
    return "undecided", sol, sol.t[-1]


def solve_m2(pot, x_lo, x_hi, n_scan=41, eps_f=1e-9, rtol=1e-12, atol=1e-15, span=200.0, n_quad=2 ** 17,
             eps_cut=1e-7):
    t_start = time.perf_counter()
    grid = np.linspace(x_lo, x_hi, n_scan)
    kinds = []
    for x in grid:
        l, d = pot.ev(x)
        kinds.append("skip" if l + 0.25 * d >= 0 else _m2_shot(pot, x, eps_f, rtol, atol, span)[0])
    br = [(grid[i], grid[i + 1]) for i in range(n_scan - 1) if kinds[i] == "under" and kinds[i + 1] in ("over", "escape")]
    if not br:
        return BounceResult("M2", False, seconds=time.perf_counter() - t_start,
                            info=dict(reason="no undershoot/overshoot bracket", scan="".join(k[0] for k in kinds)))
    a, b = br[0]
    nshot = n_scan
    for _ in range(200):
        m = 0.5 * (a + b)
        if m <= a or m >= b:
            break
        k = _m2_shot(pot, m, eps_f, rtol, atol, span)[0]
        nshot += 1
        if k == "under":
            a = m
        elif k in ("over", "escape"):
            b = m
        else:
            break
    kind, sol, t_end = _m2_shot(pot, a, eps_f, rtol, atol, span)
    tau = np.linspace(sol.t[0], t_end, 200001)
    f = sol.sol(tau)[0]
    ipk = int(np.argmax(f))
    after = np.nonzero(f[ipk:] < eps_cut * f[ipk])[0]
    t_cut = tau[ipk + after[0]] if after.size else t_end
    tq = np.linspace(sol.t[0], t_cut, n_quad + 1)
    fq, fpq = sol.sol(tq)
    w = np.full(n_quad + 1, 2.0); w[1::2] = 4.0; w[0] = w[-1] = 1.0; w *= (tq[1] - tq[0]) / 3.0
    out = _functionals(pot, tq, fq, fpq, w)
    return BounceResult("M2", True, S=out["S"], T=out["T"], V=out["V"], derrick=out["derrick"],
                        dilatation=out["dilatation"], lnphi0=a, lnR=out["lnR"],
                        seconds=time.perf_counter() - t_start,
                        info=dict(bracket_width=b - a, shots=nshot, n_brackets=len(br), f_max=out["f_max"],
                                  scan="".join(k[0] for k in kinds)))


# ----------------------------------------------------------------------------- M1: (s, u) shooting
def _m1_shot(pot, x0, rtol, atol, s_max):
    l0, d0 = pot.ev(x0)
    G0 = l0 + 0.25 * d0
    s0 = 1e-4 / math.sqrt(abs(G0))
    y0 = [1.0 + G0 * s0 ** 2 / 8.0, G0 * s0 / 4.0]
    ev = pot.ev

    def rhs(s, y):
        u, up = y
        if u <= 0.0:
            return [up, -3.0 * up / s]
        l, d = ev(x0 + math.log(u))
        return [up, -3.0 * up / s + (l + 0.25 * d) * u ** 3]

    def over(s, y):
        return y[0]
    over.terminal, over.direction = True, -1

    def under(s, y):
        return y[1]
    under.terminal, under.direction = True, 1

    def escape(s, y):
        return y[0] - 1.5
    escape.terminal, escape.direction = True, 1

    sol = solve_ivp(rhs, (s0, s_max), y0, method="DOP853", rtol=rtol, atol=atol, dense_output=True,
                    events=[over, under, escape])
    if sol.t_events[2].size:
        return "escape", sol, sol.t_events[2][0]
    if sol.t_events[0].size:
        return "over", sol, sol.t_events[0][0]
    if sol.t_events[1].size:
        return "under", sol, sol.t_events[1][0]
    return "undecided", sol, sol.t[-1]


def solve_m1(pot, x_lo, x_hi, n_scan=41, rtol=1e-12, atol=1e-18, s_max=1e16, n_quad=2 ** 17, u_cut=1e-7):
    t_start = time.perf_counter()
    grid = np.linspace(x_lo, x_hi, n_scan)
    kinds = []
    for x in grid:
        l, d = pot.ev(x)
        kinds.append("skip" if l + 0.25 * d >= 0 else _m1_shot(pot, x, rtol, atol, s_max)[0])
    br = [(grid[i], grid[i + 1]) for i in range(n_scan - 1) if kinds[i] == "under" and kinds[i + 1] in ("over", "escape")]
    if not br:
        return BounceResult("M1", False, seconds=time.perf_counter() - t_start,
                            info=dict(reason="no undershoot/overshoot bracket", scan="".join(k[0] for k in kinds)))
    a, b = br[0]
    nshot = n_scan
    for _ in range(200):
        m = 0.5 * (a + b)
        if m <= a or m >= b:
            break
        k = _m1_shot(pot, m, rtol, atol, s_max)[0]
        nshot += 1
        if k == "under":
            a = m
        elif k in ("over", "escape"):
            b = m
        else:
            break
    kind, sol, s_end = _m1_shot(pot, a, rtol, atol, s_max)
    # truncate where u has fallen to u_cut (well before the near-bounce trajectory departs)
    sg = np.geomspace(sol.t[0], s_end, 200001)
    u = sol.sol(sg)[0]
    idx = np.nonzero(u < u_cut)[0]
    s_cut = sg[idx[0]] if idx.size else s_end
    sig = np.linspace(math.log(sol.t[0]), math.log(s_cut), n_quad + 1)
    s = np.exp(sig)
    uu, up = sol.sol(s)
    x = a + np.log(np.maximum(uu, 1e-300))
    lam, beta = pot.lam(x), pot.dlam(x)
    w = np.full(n_quad + 1, 2.0); w[1::2] = 4.0; w[0] = w[-1] = 1.0; w *= (sig[1] - sig[0]) / 3.0
    C2 = (uu[-1] * s_cut ** 2) ** 2
    T = TWO_PI2 * (np.dot(w, s ** 4 * 0.5 * up ** 2) + C2 / s_cut ** 2)
    V = TWO_PI2 * np.dot(w, s ** 4 * 0.25 * lam * uu ** 4)
    dil = np.dot(w, s ** 4 * beta * uu ** 4) / np.dot(w, s ** 4 * np.abs(lam) * uu ** 4)
    lnR = float(math.log(s[int(np.argmax(s * uu))]) - a)
    return BounceResult("M1", True, S=T + V, T=T, V=V, derrick=(T + 2 * V) / T, dilatation=dil, lnphi0=a, lnR=lnR,
                        seconds=time.perf_counter() - t_start,
                        info=dict(bracket_width=b - a, shots=nshot, s_cut=s_cut, scan="".join(k[0] for k in kinds)))


# ----------------------------------------------------------------------------- M3: solve_bvp collocation
def _fl_guess(lam_ref, lnh):
    A = math.sqrt(2.0 / abs(lam_ref))
    return A, math.log(2.0 * A) - lnh          # f = A sech(tau - tau_c), phi(0) = 2 A e^{-tau_c}


def solve_m3(pot, lam_ref, lnh_guess, L=15.0, tol=1e-10, n_init=401, max_nodes=500000):
    t_start = time.perf_counter()
    A, tc = _fl_guess(lam_ref, lnh_guess)
    tau = np.linspace(tc - L, tc + L, n_init)
    sech = 1.0 / np.cosh(tau - tc)
    y = np.vstack([A * sech, -A * sech * np.tanh(tau - tc)])

    def fun(t, y):
        f = y[0]
        x = np.log(np.maximum(f, 1e-300)) - t
        G = pot.lam(x) + 0.25 * pot.dlam(x)
        return np.vstack([y[1], f + f ** 3 * G])

    def fun_jac(t, y):
        f = y[0]
        x = np.log(np.maximum(f, 1e-300)) - t
        G = pot.lam(x) + 0.25 * pot.dlam(x)
        Gx = pot.dlam(x) + 0.25 * pot.d2lam(x)
        J = np.zeros((2, 2, t.size))
        J[0, 1] = 1.0
        J[1, 0] = 1.0 + 3.0 * f ** 2 * G + f ** 2 * Gx
        return J

    def bc(ya, yb):
        return np.array([ya[1] - ya[0], yb[1] + yb[0]])

    def bc_jac(ya, yb):
        return np.array([[-1.0, 1.0], [0.0, 0.0]]), np.array([[0.0, 0.0], [1.0, 1.0]])

    sol = solve_bvp(fun, bc, tau, y, fun_jac=fun_jac, bc_jac=bc_jac, tol=tol, max_nodes=max_nodes)
    if sol.status != 0 or np.max(sol.sol(tau)[0]) < 0.1 * A:
        return BounceResult("M3", False, seconds=time.perf_counter() - t_start,
                            info=dict(reason=sol.message, nodes=sol.x.size))
    n = 2 ** 16
    tq = np.linspace(sol.x[0], sol.x[-1], n + 1)
    fq, fpq = sol.sol(tq)
    w = np.full(n + 1, 2.0); w[1::2] = 4.0; w[0] = w[-1] = 1.0; w *= (tq[1] - tq[0]) / 3.0
    out = _functionals(pot, tq, fq, fpq, w)
    x0 = float(np.log(fq[0]) - tq[0])
    return BounceResult("M3", True, S=out["S"], T=out["T"], V=out["V"], derrick=out["derrick"],
                        dilatation=out["dilatation"], lnphi0=x0, lnR=out["lnR"], seconds=time.perf_counter() - t_start,
                        info=dict(nodes=sol.x.size, rms_residual=float(np.max(sol.rms_residuals))))


# ----------------------------------------------------------------------------- M4: Chebyshev collocation
def _cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(c, 1.0 / c) / (X - X.T + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


def _clencurt(N):
    theta = np.pi * np.arange(N + 1) / N
    w = np.zeros(N + 1)
    v = np.ones(N - 1)
    if N % 2 == 0:
        w[0] = w[N] = 1.0 / (N ** 2 - 1)
        for k in range(1, N // 2):
            v -= 2.0 * np.cos(2 * k * theta[1:-1]) / (4 * k * k - 1)
        v -= np.cos(N * theta[1:-1]) / (N ** 2 - 1)
    else:
        w[0] = w[N] = 1.0 / N ** 2
        for k in range(1, (N - 1) // 2 + 1):
            v -= 2.0 * np.cos(2 * k * theta[1:-1]) / (4 * k * k - 1)
    w[1:-1] = 2.0 * v / N
    return w


def solve_m4(pot, lam_ref, lnh_guess, N=256, L=15.0, tol=1e-14, max_iter=60, init=None):
    t_start = time.perf_counter()
    D, xi = _cheb(N)
    A, tc = _fl_guess(lam_ref, lnh_guess)
    if init is not None:              # (tau_c, callable f(tau)) from another solver
        tc, finit = init
    tau = tc + L * xi
    D1 = D / L
    D2 = D1 @ D1
    f = A / np.cosh(tau - tc) if init is None else finit(tau)
    I = np.eye(N + 1)

    def residual(f):
        x = np.log(np.maximum(f, 1e-300)) - tau
        G = pot.lam(x) + 0.25 * pot.dlam(x)
        R = D2 @ f - f - f ** 3 * G
        R[0] = (D1 @ f)[0] + f[0]          # tau = tc + L (right end): f' = -f
        R[-1] = (D1 @ f)[-1] - f[-1]       # tau = tc - L (left end):  f' = +f
        return R, x, G

    R, x, G = residual(f)
    it = 0
    for it in range(1, max_iter + 1):
        Gx = pot.dlam(x) + 0.25 * pot.d2lam(x)
        J = D2 - np.diag(1.0 + 3.0 * f ** 2 * G + f ** 2 * Gx)
        J[0] = D1[0] + I[0]
        J[-1] = D1[-1] - I[-1]
        delta = np.linalg.solve(J, -R)
        step, n0 = 1.0, np.linalg.norm(R)
        while True:
            fn = f + step * delta
            Rn, xn, Gn = residual(fn)
            if np.linalg.norm(Rn) < n0 or step < 1e-4:
                break
            step *= 0.5
        f, R, x, G = fn, Rn, xn, Gn
        if np.max(np.abs(step * delta)) < tol * np.max(np.abs(f)):
            break
    ok = np.max(f) > 0.1 * A and np.all(f > 0)
    wq = _clencurt(N) * L
    fp = D1 @ f
    out = _functionals(pot, tau[::-1], f[::-1], fp[::-1], wq[::-1]) if ok else {}
    return BounceResult("M4", bool(ok), S=out.get("S", float("nan")), T=out.get("T", float("nan")),
                        V=out.get("V", float("nan")), derrick=out.get("derrick", float("nan")),
                        dilatation=out.get("dilatation", float("nan")),
                        lnphi0=float(np.log(f[-1]) - tau[-1]) if ok else float("nan"), lnR=out.get("lnR", float("nan")),
                        seconds=time.perf_counter() - t_start,
                        info=dict(N=N, L=L, newton_iterations=it, residual_inf=float(np.max(np.abs(R)))),
                        profile=(tau[::-1].copy(), f[::-1].copy()) if ok else None)


# ----------------------------------------------------------------------------- M0: legacy construction
def solve_m0(pot, x_center, half_width=8.0, n_tau=4001):
    """Legacy 'S_exact' construction applied to the same smooth potential (no phi^6, no interpolation noise):
    S(R) = 8pi^2/(3|lam_c|) + (pi^2/2) int (lambda(x(tau)) - lam_c) f_FL^4 dtau,  lam_c = lambda(x = -ln R),
    minimised over ln R. This is a first-order perturbative estimate, not a solution of the bounce equation."""
    t_start = time.perf_counter()
    s = np.linspace(-30.0, 30.0, n_tau)
    h = s[1] - s[0]
    w = np.full(n_tau, 2.0); w[1::2] = 4.0; w[0] = w[-1] = 1.0; w *= h / 3.0

    def S_of(lnR):
        lam_c = float(pot.lam(-lnR))
        if lam_c >= 0:
            return 1e100
        A = math.sqrt(2.0 / abs(lam_c))
        f = A / np.cosh(s)
        x = np.log(f) - (s + lnR)
        return 8 * math.pi ** 2 / (3 * abs(lam_c)) + 0.5 * math.pi ** 2 * np.dot(w, (pot.lam(x) - lam_c) * f ** 4)

    r = minimize_scalar(S_of, bounds=(-x_center - half_width, -x_center + half_width), method="bounded",
                        options=dict(xatol=1e-10))
    return BounceResult("M0", bool(r.fun < 1e99), S=float(r.fun), lnR=float(r.x), seconds=time.perf_counter() - t_start,
                        info=dict(note="first-order perturbative estimate; not a bounce solution"))
