import unittest

import numpy as np

from deblend_sofia_detections.support.marker_selection import (
    catalogue_marker_base,
)


class ManualMarkerSelectionTests(unittest.TestCase):
    def test_combined_mode_preserves_automatic_markers(self):
        detected = np.ma.masked_array(
            [[0, 1], [2, 0]],
            mask=[[True, False], [False, True]],
        )

        selected = catalogue_marker_base(
            detected,
            manual_markers_only=False,
        )

        self.assertIs(selected, detected)

    def test_manual_only_mode_discards_automatic_marker_labels(self):
        detected = np.ma.masked_array(
            [[0, 1], [2, 0]],
            mask=[[True, False], [False, True]],
        )

        selected = catalogue_marker_base(
            detected,
            manual_markers_only=True,
        )

        np.testing.assert_array_equal(selected.data, np.zeros((2, 2)))
        np.testing.assert_array_equal(selected.mask, detected.mask)
        self.assertIsNot(selected, detected)


if __name__ == "__main__":
    unittest.main()
