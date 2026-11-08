"""Sanity checks of audit/independent.py against analytic limits."""
import math, numpy as np
from audit import independent as I

class ConstPot:  # lambda = const < 0  -> FL: S = 8 pi^2/(3|lambda|) for any phi(0)
    def __init__(s, lam): s.l = lam; s.t_lo = 0; s.t_hi = 100; s.include_phi6 = False
    def lam(s, x): return np.full_like(np.asarray(x, float), s.l)
    def beta(s, x): return np.zeros_like(np.asarray(x, float))

class ParabolaPot:  # lambda = l0 + b (ln phi - x0)^2 : smooth minimum, exact bounce exists
    def __init__(s, l0, b, x0): s.l0, s.b, s.x0 = l0, b, x0; s.t_lo = 0; s.t_hi = 200; s.include_phi6 = False
    def lam(s, x): return s.l0 + s.b * (np.asarray(x, float) - s.x0) ** 2
    def beta(s, x): return 2 * s.b * (np.asarray(x, float) - s.x0)

p = ConstPot(-0.0123)
lb = I.lambda_bar(p, 40.0)
print("const: lambdabar =", lb, " (expect -0.0123)")
sh = I._shoot(p, 40.0, 1e-7, 1e-12, 1e-15)
act = I._action_along(p, sh.sol, sh.sol.t[0], sh.tau_end)
print("const: shot kind", sh.kind, " S =", act['S'], " exact", 8*math.pi**2/3/0.0123, " rel err", act['S']/(8*math.pi**2/3/0.0123)-1, " derrick", act['derrick'])

q = ParabolaPot(-0.0125, 2e-5, 39.0)
lo, hi = 8*math.pi**2/3/0.0125, None
lbs = [I.lambda_bar(q, x) for x in np.linspace(36, 44, 81)]
lbmin = min(lbs); up = 8*math.pi**2/3/abs(lbmin)
br, kinds = I.find_shooting_bracket(q, 30.0, 50.0, n=41)
print("parabola: brackets", br)
res = I.shoot_bounce(q, br[0][0], br[0][1])
print(f"parabola: S_lower={lo:.6f}  S_shoot={res['S']:.6f}  S_upper={up:.6f}  derrick={res['derrick']:.2e} dilatation={res['dilatation']:.2e} lnphi0={res['lnphi0']:.4f} lnphi_peak={res['lnphi_at_peak']:.4f}")
for rt in (1e-9, 1e-10, 1e-11, 1e-12):
    r2 = I.shoot_bounce(q, br[0][0], br[0][1], rtol=rt)
    print(f"   rtol={rt:g}: S={r2['S']:.9f} derrick={r2['derrick']:.1e}")
