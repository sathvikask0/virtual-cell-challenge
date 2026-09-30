import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np

from prediction_io import write_prediction_blocks


class PredictionIOTests(unittest.TestCase):
    def test_streamed_matrix_and_labels_match(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pred.h5ad"
            blocks = [np.array([[1, 0], [0, 2]]), np.array([[3, 4]])]
            write_prediction_blocks(path, ["A", "B"], ["control", "control", "A"], iter(blocks))
            a = ad.read_h5ad(path)
            np.testing.assert_array_equal(a.X.toarray(), np.vstack(blocks))
            self.assertEqual(a.obs.target.tolist(), ["control", "control", "A"])
            self.assertEqual(a.var_names.tolist(), ["A", "B"])
            self.assertTrue(a.obs_names.is_unique)

    def test_missing_rows_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                write_prediction_blocks(Path(directory) / "pred.h5ad", ["A"], ["A", "A"],
                                        iter([np.array([[1]])]))


if __name__ == "__main__":
    unittest.main()
