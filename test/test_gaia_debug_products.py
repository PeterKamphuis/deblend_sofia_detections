from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from deblend_sofia_detections.deblending.deblending import (
    write_gaia_debug_products,
)


class GaiaDebugProductTests(unittest.TestCase):
    def test_source_debug_directory_receives_masked_image_and_binary_mask(self):
        with TemporaryDirectory() as temporary_directory:
            cfg = SimpleNamespace(
                general=SimpleNamespace(debug=True),
                directories=SimpleNamespace(
                    watershed_directory=temporary_directory
                ),
            )
            wcs = WCS(naxis=2)
            wcs.wcs.crpix = [2.0, 2.0]
            wcs.wcs.cdelt = [-0.001, 0.001]
            wcs.wcs.crval = [10.0, -20.0]
            wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            gaia_mask = np.array(
                [
                    [False, False, False],
                    [False, True, True],
                    [False, False, False],
                ]
            )
            masked_background = np.arange(9, dtype=float).reshape(3, 3)
            masked_background[gaia_mask] = np.nan

            products = write_gaia_debug_products(
                cfg,
                "36",
                masked_background,
                gaia_mask,
                wcs,
                query_succeeded=True,
            )

            masked_name = Path(products["masked_background"])
            mask_name = Path(products["binary_mask"])
            self.assertEqual(
                masked_name.parent,
                Path(temporary_directory)
                / "watershed_source_36"
                / "debug_products",
            )
            self.assertTrue(masked_name.is_file())
            self.assertTrue(mask_name.is_file())
            np.testing.assert_array_equal(
                fits.getdata(mask_name), gaia_mask.astype(np.uint8)
            )
            self.assertEqual(fits.getheader(mask_name)["GAIA_OK"], True)
            self.assertEqual(fits.getheader(mask_name)["MASKNPIX"], 2)
            self.assertAlmostEqual(
                fits.getheader(mask_name)["MASKFRAC"], 2.0 / 9.0
            )
            self.assertTrue(np.isnan(fits.getdata(masked_name)[1, 1]))

    def test_debug_disabled_writes_nothing(self):
        with TemporaryDirectory() as temporary_directory:
            cfg = SimpleNamespace(
                general=SimpleNamespace(debug=False),
                directories=SimpleNamespace(
                    watershed_directory=temporary_directory
                ),
            )

            products = write_gaia_debug_products(
                cfg,
                "36",
                np.zeros((2, 2)),
                np.zeros((2, 2), dtype=bool),
                WCS(naxis=2),
            )

            self.assertEqual(products, {})
            self.assertEqual(list(Path(temporary_directory).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
