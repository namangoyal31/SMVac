"""Tests of validation/compare_approximations.py (post-processing of the committed datasets).

1. Coefficients: c1 = 13/72 is Espinosa (2020) eq. (3.21) in SMVac notation, (13 pi^2/27)/(8 pi^2/3); c2 from eq. (6.5).
2. Error statistics on synthetic data (signs, infinities, location of the maximum, identical selections).
3. Reproducibility: recomputation from the committed datasets equals the committed JSON output.
4. Thick-wall consistency: on the stored parabolic-model data the first-order error scales as c2/q^2.
"""
import json
import math
import os
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "validation"))
import compare_approximations as CA  # noqa: E402


class TestCoefficients(unittest.TestCase):
    def test_first_order_mapping(self):
        self.assertAlmostEqual(CA.C1, (13 * math.pi ** 2 / 27) / (8 * math.pi ** 2 / 3), places=15)

    def test_second_order_coefficient(self):
        # K_0 = 0.118718 is quoted in Espinosa (2020) eq. (6.6); the resulting B(5/6)/24 = -0.0136442
        self.assertAlmostEqual(CA.espinosa_second_order_coefficient(), -0.0136442, places=6)


class TestStatistics(unittest.TestCase):
    def test_signed_abs_and_infinite(self):
        snum = np.array([100.0, 200.0, 400.0, 800.0])
        est = np.array([101.0, 190.0, np.inf, 800.0])
        q = np.array([50.0, 5.0, 0.5, 0.05])
        s = CA.stats(est, snum, q, {"Mh": np.array([1.0, 2.0, 3.0, 4.0])})
        self.assertEqual(s["n"], 4)
        self.assertEqual(s["n_nonfinite"], 1)
        self.assertTrue(math.isinf(s["max_abs"]))
        self.assertAlmostEqual(s["max_abs_finite"], 0.05)
        self.assertAlmostEqual(s["signed_at_max"], -0.05)
        self.assertEqual(s["q_at_max"], 5.0)
        self.assertEqual(s["Mh"], 2.0)
        self.assertEqual(s["q_max_nonfinite"], 0.5)
        self.assertAlmostEqual(s["signed_min"], -0.05)
        self.assertAlmostEqual(s["signed_max"], 0.01)

    def test_identical_selection_for_all_estimates(self):
        df = pd.DataFrame({"S_num": [10.0, 20.0, 30.0], "S_FL": [9.0, 18.0, 25.0], "S_curv": [10.1, np.inf, 31.0],
                           "S_CGM": [10.0, 21.0, np.nan], "lambda_min": [-0.01, -0.01, -0.01],
                           "lambda2": [1e-4, 1e-3, 5e-2], "Mh": [1.0, 2.0, 3.0], "Mt": [4.0, 5.0, 6.0]})
        res = CA.compare(df, -0.0136, ["Mh", "Mt"])
        ns = {k: v["n"] for k, v in res["all"].items()}
        self.assertEqual(set(ns.values()), {3})
        self.assertEqual(res["all"]["CGM"]["n_nonfinite"], 1)      # a missing estimate counts as infinite


class TestReproducibility(unittest.TestCase):
    def test_committed_output_reproduced(self):
        path = os.path.join(ROOT, "validation", "results", "approximation_comparison.json")
        with open(path, encoding="utf-8") as fh:
            committed = json.load(fh)
        fresh = json.loads(json.dumps(CA.run(), default=float))
        self.assertEqual(fresh, committed)


class TestThickWall(unittest.TestCase):
    def test_first_order_error_scales_as_c2_over_q2(self):
        m = pd.read_csv(os.path.join(ROOT, "validation", "results", "near_critical_model.csv"))
        c2 = CA.espinosa_second_order_coefficient()
        for q in (1000.0, 300.0):
            r = m[m.q == q].iloc[0]
            e1 = r.S_FL * (1 + CA.C1 / q) / r.S_num - 1
            # S/S_FL = 1 + c1/q + c2/q^2 + O(q^-3)  =>  e1 = -c2/q^2 + O(q^-3)
            self.assertAlmostEqual(e1 * q ** 2 / -c2, 1.0, delta=5e-3 * 1000 / q)


if __name__ == "__main__":
    unittest.main()
