"""Fast tests of the independent Python reference implementation (python/smvac/reference).

Run: python -m unittest discover -s python/tests -v
"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from smvac.reference import analytic as A, bounce as B, potential as P, sm  # noqa: E402

EIGHT_PI2_3 = 8 * math.pi ** 2 / 3


class TestPhysicsInputs(unittest.TestCase):
    def test_matching_central_values(self):
        c = sm.matching(sm.Inputs(125.15, 173.34))
        self.assertAlmostEqual(c[6], 0.12604, places=12)                       # eq. (55)
        self.assertAlmostEqual(c[3], 0.93690, places=12)                       # eq. (57)
        self.assertAlmostEqual(c[1], 0.64779, places=12)                       # eq. (58)
        self.assertAlmostEqual(c[0], math.sqrt(5 / 3) * 0.35830, places=12)   # eq. (59), GUT normalisation
        self.assertAlmostEqual(c[2], 1.1666, places=12)                        # eq. (60)

    def test_power_counting(self):
        """Every beta function scales by s when couplings^2 -> s couplings^2 and (16 pi^2) -> s (16 pi^2)."""
        B2 = sm._B
        c = (0.55, 0.61, 0.83, 0.71, 0.13, 0.09, -0.031)
        saved = (B2.LOOP1, B2.LOOP2, B2.LOOP3, B2.LOOP4)
        try:
            for s in (2.0, 0.4):
                b1 = sm.beta_squares(c)
                B2.LOOP1, B2.LOOP2, B2.LOOP3, B2.LOOP4 = [saved[0] ** k * s ** k for k in (1, 2, 3, 4)]
                r = math.sqrt(s)
                b2 = sm.beta_squares((c[0] * r, c[1] * r, c[2] * r, c[3] * r, c[4] * r, c[5] * r, c[6] * s))
                B2.LOOP1, B2.LOOP2, B2.LOOP3, B2.LOOP4 = saved
                for x, y in zip(b1, b2):
                    self.assertAlmostEqual(y / (s * x), 1.0, places=12)
        finally:
            B2.LOOP1, B2.LOOP2, B2.LOOP3, B2.LOOP4 = saved

    def test_planck_scale_couplings_eq61(self):
        rg = P.RGSolution(sm.Inputs(125.15, 173.34))
        g1, g2, g3, yt, yb, ytau, lam = rg.c(rg.t_end)
        for v, ref in ((g1, 0.6154), (g2, 0.5055), (g3, 0.4873), (yt, 0.3825), (lam, -0.0143)):
            self.assertLess(abs(v - ref), 1e-4)


class TestAnalytic(unittest.TestCase):
    def test_constant(self):
        a = A.estimates(P.constant_model(-0.0125), (25, 55))
        self.assertAlmostEqual(a.S_FL / (EIGHT_PI2_3 / 0.0125), 1.0, places=14)
        self.assertAlmostEqual(a.S_CGM / (EIGHT_PI2_3 / 0.0125), 1.0, places=12)

    def test_parabola_closed_form(self):
        a = A.estimates(P.parabola_model(-0.0125, 5e-5, 39.0), (25, 55))
        closed = EIGHT_PI2_3 / (0.0125 - 13 / 72 * 5e-5)
        self.assertAlmostEqual(a.S_curv / closed, 1.0, places=13)
        self.assertAlmostEqual(a.S_CGM / closed, 1.0, places=11)

    def test_two_parameter_family(self):
        """n = 1 member = S_CGM; n -> infinity member = 4 pi^2/(l2 (q - 3/32)); finite iff q > 3/32 (model)."""
        l2 = 5e-5
        pot = P.parabola_model(-0.5 * l2, l2, 39.0)                    # q = 0.5
        a = A.estimates(pot, (15, 65))
        T, V = A.cgm2_trial(pot, a.lnh_CGM, 1.0)
        self.assertAlmostEqual(T * T / (4 * abs(V)) / a.S_CGM, 1.0, places=12)
        T, V = A.cgm2_trial(pot, 39.25, 1e5)
        self.assertAlmostEqual(T * T / (4 * abs(V)) / (4 * math.pi ** 2 / (l2 * (0.5 - 3 / 32))), 1.0, places=4)
        S2, _, n2 = A.cgm2_estimate(pot, a.lnh_CGM)
        self.assertLess(S2, a.S_CGM)
        self.assertGreater(n2, 1.0)
        for q, finite in ((0.1, True), (0.085, False)):               # thresholds 13/72 (CGM) and 3/32 (CGM2)
            pot = P.parabola_model(-q * l2, l2, 39.0)
            a = A.estimates(pot, (15, 65))
            self.assertFalse(math.isfinite(a.S_CGM))
            S2, _, _ = A.cgm2_estimate(pot, a.x_min + 5.0 / 6.0)
            self.assertEqual(math.isfinite(S2), finite)


class TestBounce(unittest.TestCase):
    def test_spectral_vs_reference(self):
        pot = P.parabola_model(-0.0125, 5e-5, 39.0)
        a = A.estimates(pot, (25, 55))
        r = B.solve_m4(pot, a.lambda_min, a.lnh_CGM, N=256)
        self.assertTrue(r.ok)
        self.assertAlmostEqual(r.S / 2107.0357973921, 1.0, places=10)   # M1/M2/M3/M4 consensus
        self.assertLess(abs(r.derrick), 1e-10)
        self.assertLess(abs(r.dilatation), 1e-10)
        self.assertLessEqual(a.S_FL, r.S)
        self.assertLessEqual(r.S, a.S_CGM)

    def test_shooting_agrees_with_spectral_thin_wall(self):
        pot = P.parabola_model(-2.5 * 5e-5, 5e-5, 39.0)
        a = A.estimates(pot, (25, 55))
        lo, hi = B.admissible_interval(pot, a.x_min)
        r2 = B.solve_m2(pot, lo, hi, n_scan=17)
        r4 = B.solve_m4(pot, a.lambda_min, a.lnh_CGM, N=384)
        self.assertTrue(r2.ok and r4.ok)
        self.assertAlmostEqual(r2.S / r4.S, 1.0, places=9)

    def test_legacy_construction_unbounded(self):
        """The legacy S_exact construction is unbounded below even on a smooth potential (q = 25)."""
        pot = P.parabola_model(-25 * 5e-5, 5e-5, 39.0)
        a = A.estimates(pot, (25, 55))
        self.assertLess(B.solve_m0(pot, a.x_min).S, 0.0)


if __name__ == "__main__":
    unittest.main()
