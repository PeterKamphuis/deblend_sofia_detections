import importlib.util
from pathlib import Path
import unittest

import numpy as np


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "convert_optical_fits_to_2d.py"
)
SPEC = importlib.util.spec_from_file_location(
    "convert_optical_fits_to_2d", SCRIPT_PATH
)
CONVERTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONVERTER)


class StandaloneOpticalConverterTests(unittest.TestCase):
    def test_default_output_name_is_created_beside_input(self):
        output = CONVERTER.default_output_path("/data/field.fits")

        self.assertEqual(output, Path("/data/field_2d.fits"))

    def test_default_output_name_preserves_compressed_fits_suffix(self):
        output = CONVERTER.default_output_path("/data/field.fits.gz")

        self.assertEqual(output, Path("/data/field_2d.fits.gz"))

    def test_explicit_output_argument_remains_optional_in_parser(self):
        arguments = CONVERTER.build_argument_parser().parse_args(
            ["input.fits"]
        )

        self.assertEqual(arguments.input_fits, "input.fits")
        self.assertIsNone(arguments.output_fits)

    def test_explicit_output_argument_is_still_accepted(self):
        arguments = CONVERTER.build_argument_parser().parse_args(
            ["input.fits", "custom.fits"]
        )

        self.assertEqual(arguments.output_fits, "custom.fits")

    def test_arbitrary_leading_dimensions_are_collapsed(self):
        image = np.arange(2 * 3 * 4 * 5, dtype=float).reshape(2, 3, 4, 5)

        reduced, axes = CONVERTER.collapse_to_2d(
            image, spatial_axes=(2, 3), method="mean"
        )

        np.testing.assert_allclose(reduced, np.mean(image, axis=(0, 1)))
        self.assertEqual(axes, (0, 1))

    def test_first_method_selects_first_non_spatial_plane(self):
        image = np.arange(3 * 2 * 4, dtype=float).reshape(3, 2, 4)

        reduced, axes = CONVERTER.collapse_to_2d(
            image, spatial_axes=(1, 2), method="first"
        )

        np.testing.assert_array_equal(reduced, image[0])
        self.assertEqual(axes, (0,))

    def test_rgb_axis_is_inferred_when_wcs_mapping_is_unavailable(self):
        class IncompleteWCS:
            pixel_n_dim = 2

            def get_axis_types(self):
                return []

        image = np.ones((2, 5, 3))

        axes, source = CONVERTER.infer_spatial_axes(image, IncompleteWCS())

        self.assertEqual(axes, (0, 1))
        self.assertEqual(source, "RGB(A) axis length")

    def test_input_cannot_be_reduced_from_one_dimension(self):
        with self.assertRaisesRegex(ValueError, "at least two dimensions"):
            CONVERTER.collapse_to_2d(
                np.ones(5), spatial_axes=(0, 0), method="mean"
            )


if __name__ == "__main__":
    unittest.main()
