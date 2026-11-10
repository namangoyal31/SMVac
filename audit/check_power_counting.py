"""Loop-order power-counting check of the legacy beta functions (literature-independent).

Assign degree 1 to each of g_i^2, y_i^2, lambda.  An n-loop contribution to d(x^2)/d ln mu^2
(x = gauge or Yukawa coupling) or to d(lambda)/d ln mu^2 must have total degree n+1 and carry
(16 pi^2)^-n.  Hence, scaling couplings^2 -> s * couplings^2 together with (16pi^2)^n -> s^n (16pi^2)^n
must scale a correct beta function exactly by s.  We test each loop order separately by isolating
the 1/LOOPn pieces (setting the other LOOP factors to infinity)."""
import math, random
from audit import _legacy_betas_extracted as B

NAMES = ['betaG1sq', 'betaG2sq', 'betaG3sq', 'betaLambda', 'betaYt2', 'betaYb2', 'betaYtau2']
BIG = 1e300

def term(fn, c, n, s):
    """n-loop piece of beta at couplings^2 scaled by s, with LOOPn scaled by s^n."""
    g1, g2, g3, yt, yb, ytau, lam = c
    r = math.sqrt(s)
    for k in (1, 2, 3, 4):
        setattr(B, f'LOOP{k}', (s * 157.91367041742973) ** k if k == n else BIG)
    return fn(g1 * r, g2 * r, g3 * r, yt * r, yb * r, ytau * r, lam * s)

random.seed(1)
bad = []
for name in NAMES:
    fn = getattr(B, name)
    for n in (1, 2, 3, 4):
        worst = 0.0
        for _ in range(20):
            c = [random.uniform(0.2, 1.2) for _ in range(6)] + [random.uniform(-0.3, 0.3)]
            t1, t2 = term(fn, c, n, 1.0), term(fn, c, n, 2.0)
            if abs(t1) < 1e-200 and abs(t2) < 1e-200:  # loop order absent in this function
                continue
            worst = max(worst, abs(t2 / 2.0 - t1) / max(abs(t1), 1e-300))
        status = ("absent" if worst == 0.0 else "OK") if worst < 1e-12 else "INCONSISTENT"
        if worst > 0 or n <= 3:
            print(f"{name:11s} {n}-loop: max relative violation {worst:.3e}  {status}")
        if status == "INCONSISTENT":
            bad.append((name, n))
print("INCONSISTENT pieces:", bad)
