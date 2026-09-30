import unittest

import numpy as np

from learned_transfer import evaluate, predict_correction, train_correction


class LearnedTransferTests(unittest.TestCase):
    def test_identical_line_effects_need_no_correction(self):
        genes = np.array(["A", "B", "C", "D"])
        targets = np.array(["A", "B", "C"])
        counts = np.array([[4., 8., 12., 20.], [10., 3., 15., 12.], [9., 16., 2., 13.]])
        ctrl = np.array([10., 10., 10., 10.])
        train = {name: dict(targets=targets, counts=counts * depth, ctrl=ctrl * depth,
                            tidx={t: i for i, t in enumerate(targets)})
                 for name, depth in [("first", 1.), ("second", 2.)]}
        model = train_correction(train, genes, rank=2)
        base, residual = predict_correction(model, train, targets, ctrl)
        np.testing.assert_allclose(residual, 0., atol=1e-12)
        self.assertEqual(base.shape, (3, 4))
        self.assertEqual(model["train_pairs"], 6)

    def test_control_prediction_has_chance_discrimination(self):
        genes = np.array(["A", "B", "C", "D"])
        targets = ["A", "B", "C"]
        ctrl = np.array([1., 2., 3., 4.])
        true = np.array([[.1, 2., 3., 5.], [1., .1, 3., 3.], [1., 2., .1, 6.]])
        pred = np.repeat(ctrl[None], 3, axis=0)
        scores = evaluate(pred, true, ctrl, targets, genes)
        self.assertAlmostEqual(scores["relative_logbulk_mse"], 1.)
        self.assertAlmostEqual(scores["pds"], .5)


if __name__ == "__main__":
    unittest.main()
