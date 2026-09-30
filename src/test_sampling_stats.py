import unittest

import numpy as np

from predict_2026 import context_stats, sample_cells
from sampling_stats import depth_aware_stats


class SamplingStatsTests(unittest.TestCase):
    def test_variable_depth_simulation_recovers_dispersion_better(self):
        true_phi = np.full(100, .3)
        X = sample_cells(np.full(100, .01), np.array([250., 2500.]), true_phi, 0.,
                         np.random.default_rng(42), n=10000)
        _, _, old_phi = context_stats(X)
        share, totals, corrected_phi = depth_aware_stats(X)
        self.assertLess(np.median(np.abs(corrected_phi - true_phi)),
                        np.median(np.abs(old_phi - true_phi)))
        self.assertAlmostEqual(share.sum(), 1.)
        self.assertTrue(np.all(totals > 0))

    def test_zero_depth_rejected(self):
        with self.assertRaises(ValueError):
            depth_aware_stats(np.zeros((2, 3)))


if __name__ == "__main__":
    unittest.main()
