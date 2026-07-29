import unittest

import numpy as np

from deblend_sofia_detections.support.optical_image import (
    celestial_numpy_axes,
    collapse_optical_data,
)


class CollapseOpticalDataTests(unittest.TestCase):
    def test_standard_fits_wcs_axes_map_to_y_and_x_array_axes(self):
        axis_types = [
            {"coordinate_type": "celestial"},
            {"coordinate_type": "celestial"},
            {"coordinate_type": "scalar"},
        ]
        correlation = np.eye(3, dtype=bool)

        axes = celestial_numpy_axes(axis_types, correlation, 3, 3)

        self.assertEqual(axes, (1, 2))

    def test_channel_last_layout_is_recovered_from_wcs_correlation(self):
        axis_types = [
            {"coordinate_type": "scalar"},
            {"coordinate_type": "celestial"},
            {"coordinate_type": "celestial"},
        ]
        correlation = np.eye(3, dtype=bool)

        axes = celestial_numpy_axes(axis_types, correlation, 3, 3)

        self.assertEqual(axes, (0, 1))

    def test_two_dimensional_image_is_unchanged(self):
        image = np.arange(6).reshape(2, 3)

        reduced, collapsed_axes = collapse_optical_data(image)

        np.testing.assert_array_equal(reduced, image)
        self.assertEqual(collapsed_axes, ())

    def test_standard_rgb_first_image_is_averaged_to_grayscale(self):
        image = np.stack(
            [
                np.full((2, 2), 1.0),
                np.full((2, 2), 4.0),
                np.full((2, 2), 7.0),
            ]
        )

        reduced, collapsed_axes = collapse_optical_data(image)

        np.testing.assert_allclose(reduced, np.full((2, 2), 4.0))
        self.assertEqual(collapsed_axes, (0,))

    def test_wcs_selected_channel_last_axis_is_collapsed(self):
        image = np.stack(
            [
                np.full((2, 3), 2.0),
                np.full((2, 3), 6.0),
                np.full((2, 3), 10.0),
            ],
            axis=-1,
        )

        reduced, collapsed_axes = collapse_optical_data(
            image, celestial_numpy_axes=(0, 1)
        )

        np.testing.assert_allclose(reduced, np.full((2, 3), 6.0))
        self.assertEqual(collapsed_axes, (2,))

    def test_singleton_and_colour_axes_can_both_be_collapsed_from_wcs(self):
        image = np.ones((1, 3, 2, 4))

        reduced, collapsed_axes = collapse_optical_data(
            image, celestial_numpy_axes=(2, 3)
        )

        self.assertEqual(reduced.shape, (2, 4))
        self.assertEqual(collapsed_axes, (0, 1))

    def test_fewer_than_two_dimensions_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "at least two dimensions"):
            collapse_optical_data(np.ones(5))


if __name__ == "__main__":
    unittest.main()
