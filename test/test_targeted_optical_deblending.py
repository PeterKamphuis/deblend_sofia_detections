from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from deblend_sofia_detections.catalogue.dr10 import (
    optical_labels_with_multiple_moment0_peaks,
)
from deblend_sofia_detections.deblending.deblending import (
    targeted_deblend_optical_regions,
    write_targeted_optical_deblend_debug_products,
)
from deblend_sofia_detections.support.errors import InputError


def make_wcs(width=64, height=64):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [width / 2.0, height / 2.0]
    wcs.wcs.cdelt = np.array([-0.001, 0.001])
    wcs.wcs.crval = [86.48, -25.75]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = wcs.to_header()
    header["NAXIS1"] = width
    header["NAXIS2"] = height
    return wcs, header


def gaussian_image(shape, components):
    y_grid, x_grid = np.indices(shape)
    image = np.zeros(shape, dtype=float)
    for x_centre, y_centre, amplitude, sigma in components:
        image += amplitude * np.exp(
            -(
                (x_grid - x_centre) ** 2
                + (y_grid - y_centre) ** 2
            )
            / (2.0 * sigma**2)
        )
    return image


class TargetedOpticalDeblendingTests(unittest.TestCase):
    def test_only_label_with_multiple_moment0_peaks_is_selected(self):
        _, header = make_wcs()
        markers = np.zeros((64, 64), dtype=int)
        markers[18:47, 10:47] = 1
        markers[48:60, 48:60] = 2
        peak_map = np.zeros((64, 64), dtype=np.uint8)
        peak_map[32, 20] = 1
        peak_map[32, 36] = 1
        peak_map[54, 54] = 1

        labels = optical_labels_with_multiple_moment0_peaks(
            peak_map, header, markers, header
        )

        self.assertEqual(labels, [1])

    def test_source93_like_region_is_split_while_other_label_is_unchanged(self):
        _, header = make_wcs()
        image = gaussian_image(
            (64, 64),
            [
                (20, 32, 100.0, 5.0),
                (36, 32, 70.0, 5.0),
                (54, 54, 40.0, 2.0),
            ],
        )
        markers = np.zeros((64, 64), dtype=int)
        markers[18:47, 10:47] = 1
        markers[48:60, 48:60] = 2
        peak_map = np.zeros((64, 64), dtype=np.uint8)
        peak_map[32, 20] = 1
        peak_map[32, 36] = 1
        peak_map[54, 54] = 1

        deblended, target_labels = targeted_deblend_optical_regions(
            image,
            markers,
            header,
            peak_map,
            header,
            nlevels=32,
            contrast=0.001,
            min_pixels=5,
        )
        result = deblended.filled(0)

        self.assertEqual(target_labels, [1])
        self.assertNotEqual(int(result[32, 20]), int(result[32, 36]))
        label_two_children = np.unique(result[markers == 2])
        self.assertEqual(list(label_two_children[label_two_children > 0]), [1])
        self.assertEqual(
            len(np.unique(result[result > 0])),
            3,
            "two blended peaks plus the untouched region are expected",
        )

    def test_source36_like_separate_peak_labels_are_preserved_exactly(self):
        _, header = make_wcs()
        image = gaussian_image(
            (64, 64),
            [(18, 30, 80.0, 4.0), (44, 30, 70.0, 4.0)],
        )
        markers = np.zeros((64, 64), dtype=int)
        markers[18:43, 7:29] = 1
        markers[18:43, 34:56] = 2
        peak_map = np.zeros((64, 64), dtype=np.uint8)
        peak_map[30, 18] = 1
        peak_map[30, 44] = 1

        result, target_labels = targeted_deblend_optical_regions(
            image,
            markers,
            header,
            peak_map,
            header,
        )

        self.assertEqual(target_labels, [])
        np.testing.assert_array_equal(result.filled(0), markers)

    def test_invalid_targeted_deblend_parameters_are_rejected(self):
        _, header = make_wcs(width=8, height=8)
        data = np.ones((8, 8), dtype=float)
        markers = np.ones((8, 8), dtype=int)
        peaks = np.zeros((8, 8), dtype=np.uint8)

        with self.assertRaisesRegex(InputError, "nlevels"):
            targeted_deblend_optical_regions(
                data,
                markers,
                header,
                peaks,
                header,
                nlevels=1,
            )
        with self.assertRaisesRegex(InputError, "contrast"):
            targeted_deblend_optical_regions(
                data,
                markers,
                header,
                peaks,
                header,
                contrast=1.1,
            )

    def test_debug_products_preserve_before_and_after_segmentations(self):
        _, header = make_wcs(width=8, height=8)
        before = np.zeros((8, 8), dtype=int)
        before[1:7, 1:7] = 1
        after = before.copy()
        after[1:7, 4:7] = 2

        with TemporaryDirectory() as temporary_directory:
            before_name, after_name = (
                write_targeted_optical_deblend_debug_products(
                    temporary_directory,
                    "93",
                    header,
                    before,
                    after,
                    [1],
                    nlevels=32,
                    contrast=0.001,
                    min_pixels=20,
                )
            )

            np.testing.assert_array_equal(fits.getdata(before_name), before)
            np.testing.assert_array_equal(fits.getdata(after_name), after)
            output_header = fits.getheader(after_name)
            self.assertTrue(output_header["TDBLEND"])
            self.assertEqual(output_header["TDBNRAW"], 1)
            self.assertEqual(output_header["TDBNFIN"], 2)
            self.assertEqual(
                Path(after_name).name,
                "optical_segmentation_after_targeted_deblend_source_93.fits",
            )


if __name__ == "__main__":
    unittest.main()
