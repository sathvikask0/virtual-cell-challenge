import unittest

import numpy as np

from correlated_noise import fit_factors, sample_lognormal


class CorrelatedNoiseTests(unittest.TestCase):
    def test_marginal_moments_and_shared_covariance(self):
        phi = np.array([.2, .3])
        loading = np.array([[np.sqrt(.1)], [np.sqrt(.15)]])
        X = sample_lognormal(np.array([.4, .6]), np.array([50.]), phi, 0.,
                             np.random.default_rng(10), n=100000, factors=loading).toarray()
        mean = np.array([20., 30.])
        np.testing.assert_allclose(X.mean(0), mean, rtol=.015)
        np.testing.assert_allclose(X.var(0), mean + phi * mean ** 2, rtol=.04)
        expected_cov = np.prod(mean) * np.expm1(np.prod(loading))
        self.assertAlmostEqual(np.cov(X.T)[0, 1] / expected_cov, 1., delta=.05)

    def test_no_variation_controls_produce_no_factors(self):
        factors = fit_factors(np.ones((4, 3)), np.zeros(3))
        self.assertEqual(factors.shape, (3, 0))

    def test_excess_shared_variance_rejected(self):
        with self.assertRaises(ValueError):
            sample_lognormal(np.array([1.]), np.array([10.]), np.array([.1]), 0.,
                             np.random.default_rng(0), factors=np.array([[1.]]))


if __name__ == "__main__":
    unittest.main()
