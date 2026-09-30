"""Behavior checks for calibration invariants; run with unittest discovery."""
import unittest

import numpy as np

from calibration import calibrated_change, mean_prediction
from predict_2026 import target_change


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.gidx = {"A": 0, "B": 1, "C": 2}
        self.share = np.array([.2, .3, .5])
        self.src = ({"A": np.array([-.8, .2, -.1])}, {"A": np.ones(3)},
                    {"A": [-.2]}, np.array([-.1, .1, 0.]), -.05, -.9)

    def test_defaults_preserve_existing_transfer(self):
        for target in ["A", "B", "unmeasured"]:
            expected = target_change(target, self.gidx, self.src, .5)
            actual = calibrated_change(target, self.gidx, self.src, self.share)
            np.testing.assert_allclose(actual[0], expected[0])
            self.assertEqual(actual[1], expected[1])

    def test_no_change_is_control_even_for_seen_target(self):
        pred = mean_prediction("A", self.gidx, self.src, self.share, 100., mode="control")
        np.testing.assert_allclose(pred, [20., 30., 50.])

    def test_total_factor_changes_depth_without_changing_composition(self):
        a = mean_prediction("A", self.gidx, self.src, self.share, 100., total_alpha=0.)
        b = mean_prediction("A", self.gidx, self.src, self.share, 100., total_alpha=1.)
        np.testing.assert_allclose(a / a.sum(), b / b.sum())
        self.assertAlmostEqual(a.sum(), 100.)
        self.assertAlmostEqual(b.sum(), 100 * np.exp(-.2))

    def test_target_only_and_disabled_target_drop(self):
        pred = mean_prediction("A", self.gidx, self.src, self.share, 100.,
                               mode="target-only", own_alpha=0.)
        np.testing.assert_allclose(pred, self.share * 100.)


if __name__ == "__main__":
    unittest.main()
