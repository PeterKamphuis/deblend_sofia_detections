import io
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch
from urllib.parse import parse_qs, urlparse

import numpy as np
from astropy.io import fits
from astropy.table import QTable
from astropy.wcs import WCS

from deblend_sofia_detections.catalogue.dr10 import (
    DR10_COLUMNS,
    build_dr10_query,
    detect_positive_beam_scale_moment0_peaks,
    download_dr10_catalogue,
    filter_dr10_counterparts_by_moment0_peaks,
    normalise_galaxy_types,
    prepare_dr10_catalogue,
    remove_rejected_optical_regions,
    select_dr10_counterparts,
    write_moment0_peak_filter_debug_products,
)
from deblend_sofia_detections.catalogue.search import search_counter_part
from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.deblending.image_manipulation import (
    mask_source_from_table,
)
from deblend_sofia_detections.deblending.deblending import set_optical_markers


def make_wcs(width=8, height=8):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [width / 2.0, height / 2.0]
    wcs.wcs.cdelt = np.array([-0.01, 0.01])
    wcs.wcs.crval = [10.0, -20.0]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = wcs.to_header()
    header["NAXIS1"] = width
    header["NAXIS2"] = height
    return wcs, header


def add_beam(header, bmaj=0.03, bmin=0.02, bpa=25.0):
    header = header.copy()
    header["BMAJ"] = bmaj
    header["BMIN"] = bmin
    header["BPA"] = bpa
    return header


def make_selected_table(wcs, candidates):
    """Create selected rows from (name, x, y, label, flux, type) tuples."""
    coordinates = [
        wcs.pixel_to_world_values(candidate[1], candidate[2])
        for candidate in candidates
    ]
    return QTable(
        {
            "Name": [candidate[0] for candidate in candidates],
            "RA": [coordinate[0] for coordinate in coordinates],
            "DEC": [coordinate[1] for coordinate in coordinates],
            "type": [candidate[5] for candidate in candidates],
            "flux_g": [candidate[4] for candidate in candidates],
            "optical_label": [candidate[3] for candidate in candidates],
        },
        units={"RA": "deg", "DEC": "deg"},
    )


def make_cfg(directory, *, manual_tables=None, original_tables=False):
    return SimpleNamespace(
        input=SimpleNamespace(
            manual_input_tables=(
                [None] if manual_tables is None else manual_tables
            ),
            auto_query_catalogue=True,
            galaxy_types=["REX", "EXP", "DEV", "SER"],
            original_tables=original_tables,
        ),
        general=SimpleNamespace(verbose=False),
        directories=SimpleNamespace(
            ancillary_directory=str(Path(directory) / "ancillary")
        ),
        sofia=SimpleNamespace(directory=str(directory), basename="field"),
        internal=SimpleNamespace(auto_catalogue_path=None),
    )


class Dr10CatalogueTests(unittest.TestCase):
    def test_query_contains_bounds_and_configured_types(self):
        query = build_dr10_query(
            {
                "ra_start": 85.0,
                "ra_end": 88.0,
                "ra_wraps": False,
                "dec_min": -27.0,
                "dec_max": -24.0,
            },
            ["rex", "SER"],
        )

        self.assertIn("FROM ls_dr10.tractor", query)
        self.assertIn("brick_primary = 1", query)
        self.assertIn("type IN ('REX', 'SER')", query)
        self.assertIn("ra BETWEEN 85.000000000000 AND 88.000000000000", query)
        self.assertIn("dec BETWEEN -27.000000000000 AND -24.000000000000", query)

    def test_invalid_or_empty_galaxy_types_are_rejected(self):
        with self.assertRaisesRegex(InputError, "at least one"):
            normalise_galaxy_types([])
        with self.assertRaisesRegex(InputError, "only letters"):
            normalise_galaxy_types(["SER'); DROP TABLE"])

    def test_manual_catalogue_prevents_automatic_download(self):
        cfg = make_cfg("/unused", manual_tables=["manual.csv"])
        with patch(
            "deblend_sofia_detections.catalogue.dr10.download_dr10_catalogue"
        ) as downloader:
            result = prepare_dr10_catalogue(cfg)

        self.assertIsNone(result)
        self.assertIsNone(cfg.internal.auto_catalogue_path)
        downloader.assert_not_called()

    def test_download_is_cached_with_query_metadata(self):
        with TemporaryDirectory() as temporary_directory:
            cfg = make_cfg(temporary_directory)
            _, header = make_wcs()
            fits.PrimaryHDU(
                np.zeros((8, 8)), header=header
            ).writeto(
                Path(temporary_directory) / "field_mom0.fits",
                output_verify="silentfix",
            )
            csv_response = (
                ",".join(DR10_COLUMNS)
                + "\n10000,1,1000m200,2,1,3,SER,10.0,-20.0,12.5\n"
            ).encode("utf-8")
            opener = Mock(return_value=io.BytesIO(csv_response))

            catalogue_path = Path(
                download_dr10_catalogue(cfg, urlopen_func=opener)
            )

            self.assertTrue(catalogue_path.is_file())
            self.assertTrue(catalogue_path.with_suffix(".json").is_file())
            request_query = parse_qs(
                urlparse(opener.call_args.args[0]).query
            )["QUERY"][0]
            self.assertIn("type IN ('REX', 'EXP', 'DEV', 'SER')", request_query)

            cached_opener = Mock(side_effect=AssertionError("unexpected query"))
            cached_path = download_dr10_catalogue(
                cfg, urlopen_func=cached_opener
            )

            self.assertEqual(Path(cached_path), catalogue_path)
            cached_opener.assert_not_called()

    def test_highest_flux_allowed_type_wins_per_cyan_region(self):
        wcs, header = make_wcs()
        markers = np.zeros((8, 8), dtype=int)
        markers[1:4, 1:4] = 1
        markers[5:7, 5:7] = 2
        hi_mask = markers > 0
        coordinates = [
            wcs.pixel_to_world_values(1.25, 1.25),
            wcs.pixel_to_world_values(2.25, 2.25),
            wcs.pixel_to_world_values(2.5, 1.5),
            wcs.pixel_to_world_values(5.25, 5.25),
        ]
        table = QTable(
            {
                "release": [10000] * 4,
                "brickid": [1] * 4,
                "brickname": ["1000m200"] * 4,
                "objid": [10, 11, 12, 13],
                "brick_primary": [1] * 4,
                "ls_id": [100, 101, 102, 103],
                "type": ["REX", "EXP", "PSF", "SER"],
                "ra": [coordinate[0] for coordinate in coordinates],
                "dec": [coordinate[1] for coordinate in coordinates],
                "flux_g": [5.0, 12.0, 1000.0, 7.0],
            }
        )
        # The downloaded catalogue loader attaches degree units before the
        # per-source selection runs.
        table["ra"].unit = "deg"
        table["dec"].unit = "deg"

        selected = select_dr10_counterparts(
            table,
            markers,
            header,
            hi_mask,
            ["REX", "EXP", "DEV", "SER"],
        )

        self.assertEqual(list(selected["objid"]), [11, 13])
        self.assertEqual(list(selected["optical_label"]), [1, 2])
        self.assertEqual(
            list(selected["Name"]),
            ["LS_1000m200_11", "LS_1000m200_13"],
        )

    def test_catalogue_row_outside_hi_is_not_selected(self):
        wcs, header = make_wcs()
        markers = np.zeros((8, 8), dtype=int)
        markers[1:4, 1:4] = 1
        hi_mask = np.zeros_like(markers, dtype=bool)
        ra, dec = wcs.pixel_to_world_values(2.25, 2.25)
        table = QTable(
            {
                "type": ["SER"],
                "ra": [ra],
                "dec": [dec],
                "flux_g": [20.0],
            }
        )

        selected = select_dr10_counterparts(
            table, markers, header, hi_mask, ["SER"]
        )

        self.assertEqual(len(selected), 0)

    def test_catalogue_position_uses_nearest_segmentation_pixel(self):
        wcs, header = make_wcs()
        markers = np.zeros((8, 8), dtype=int)
        markers[3, 3] = 1
        markers[3, 4] = 2
        hi_mask = markers > 0
        # This coordinate is centred closer to pixel x=4 than x=3. Integer
        # truncation would incorrectly associate it with optical label 1.
        ra, dec = wcs.pixel_to_world_values(3.75, 3.0)
        table = QTable(
            {
                "type": ["SER"],
                "ra": [ra],
                "dec": [dec],
                "flux_g": [20.0],
            },
            units={"ra": "deg", "dec": "deg"},
        )

        selected = select_dr10_counterparts(
            table, markers, header, hi_mask, ["SER"]
        )

        self.assertEqual(len(selected), 1)
        self.assertEqual(int(selected["optical_label"][0]), 2)

    def test_selected_dr10_row_can_seed_marker_without_size_columns(self):
        wcs, header = make_wcs()
        markers = np.ma.masked_array(
            np.zeros((8, 8), dtype=int), mask=np.zeros((8, 8), dtype=bool)
        )
        ra, dec = wcs.pixel_to_world_values(3.25, 3.25)
        selected = QTable(
            {
                "Name": ["LS_1000m200_42"],
                "RA": [ra],
                "DEC": [dec],
                "type": ["SER"],
                "flux_g": [10.0],
            }
        )
        cfg = SimpleNamespace(
            general=SimpleNamespace(verbose=False, debug=False)
        )

        result = mask_source_from_table(
            cfg,
            markers,
            header,
            mask=np.ones((8, 8), dtype=bool),
            src_table=selected,
        )

        self.assertGreater(result[3, 3], 0)

    def test_moment0_filter_retains_each_region_with_a_positive_peak(self):
        wcs, optical_header = make_wcs(width=20, height=20)
        moment0_header = add_beam(optical_header)
        markers = np.zeros((20, 20), dtype=int)
        markers[2:8, 2:8] = 1
        markers[2:8, 12:18] = 2
        markers[12:18, 2:8] = 3
        selected = make_selected_table(
            wcs,
            [
                ("LS_A", 5, 5, 1, 100.0, "SER"),
                ("LS_B", 15, 5, 2, 10.0, "EXP"),
                ("LS_C", 5, 15, 3, 1.0, "REX"),
            ],
        )
        moment0 = np.zeros((20, 20), dtype=float)
        moment0[5, 5] = 8.0
        # Deliberately faint: the filter has no amplitude or S/N threshold.
        moment0[5, 15] = 1e-12

        accepted, audit, peak_map = (
            filter_dr10_counterparts_by_moment0_peaks(
                selected,
                markers,
                optical_header,
                moment0,
                moment0_header,
                np.ones((3, 20, 20), dtype=bool),
            )
        )

        self.assertEqual(list(accepted["Name"]), ["LS_A", "LS_B"])
        self.assertEqual(
            list(audit["moment0_peak_supported"]), [True, True, False]
        )
        self.assertEqual(
            list(audit["status"]), ["accepted", "accepted", "rejected"]
        )
        self.assertEqual(int(np.sum(peak_map)), 2)
        self.assertAlmostEqual(float(audit["peak_value"][1]), 1e-12)
        self.assertEqual(audit["rejection_reason"][0], "")
        self.assertIn("no_positive", audit["rejection_reason"][2])

    def test_flat_topped_maximum_is_collapsed_to_one_peak(self):
        _, header = make_wcs(width=16, height=16)
        header = add_beam(header, bmaj=0.05, bmin=0.04, bpa=0.0)
        moment0 = np.zeros((16, 16), dtype=float)
        moment0[6:9, 6:9] = 2.5

        peak_map = detect_positive_beam_scale_moment0_peaks(
            moment0, header, np.ones_like(moment0, dtype=bool)
        )

        self.assertEqual(int(np.sum(peak_map)), 1)
        peak_y, peak_x = np.nonzero(peak_map)
        self.assertIn(int(peak_x[0]), (7,))
        self.assertIn(int(peak_y[0]), (7,))

    def test_peak_outside_exact_cyan_label_does_not_support_marker(self):
        wcs, header = make_wcs(width=16, height=16)
        header = add_beam(header)
        markers = np.zeros((16, 16), dtype=int)
        markers[3:8, 3:8] = 4
        selected = make_selected_table(
            wcs, [("LS_EDGE", 5, 5, 4, 4.0, "SER")]
        )
        moment0 = np.zeros((16, 16), dtype=float)
        # y=8 is immediately outside the exact [3:8] segmentation label.
        moment0[8, 5] = 9.0

        accepted, audit, _ = filter_dr10_counterparts_by_moment0_peaks(
            selected,
            markers,
            header,
            moment0,
            header,
            np.ones_like(moment0, dtype=bool),
        )

        self.assertEqual(len(accepted), 0)
        self.assertFalse(bool(audit["moment0_peak_supported"][0]))

    def test_peak_mapping_supports_different_wcs_resolutions(self):
        optical_wcs, optical_header = make_wcs(width=20, height=20)
        moment0_wcs, moment0_header = make_wcs(width=10, height=10)
        moment0_wcs.wcs.cdelt = np.array([-0.02, 0.02])
        moment0_header = add_beam(moment0_wcs.to_header(), 0.04, 0.03, 45.0)
        moment0_header["NAXIS1"] = 10
        moment0_header["NAXIS2"] = 10

        peak_pixel = (7, 6)
        peak_ra, peak_dec = moment0_wcs.pixel_to_world_values(*peak_pixel)
        optical_x, optical_y = optical_wcs.world_to_pixel_values(
            peak_ra, peak_dec
        )
        optical_x = int(np.rint(optical_x))
        optical_y = int(np.rint(optical_y))
        markers = np.zeros((20, 20), dtype=int)
        markers[
            optical_y - 1 : optical_y + 2,
            optical_x - 1 : optical_x + 2,
        ] = 7
        selected = make_selected_table(
            optical_wcs,
            [("LS_RES", optical_x, optical_y, 7, 3.0, "DEV")],
        )
        moment0 = np.zeros((10, 10), dtype=float)
        moment0[peak_pixel[1], peak_pixel[0]] = 1.0

        accepted, audit, _ = filter_dr10_counterparts_by_moment0_peaks(
            selected,
            markers,
            optical_header,
            moment0,
            moment0_header,
            np.ones_like(moment0, dtype=bool),
        )

        self.assertEqual(list(accepted["Name"]), ["LS_RES"])
        self.assertTrue(bool(audit["moment0_peak_supported"][0]))

    def test_invalid_moment0_wcs_or_beam_produces_input_error(self):
        _, valid_header = make_wcs(width=8, height=8)
        data = np.ones((8, 8), dtype=float)
        mask = np.ones_like(data, dtype=bool)

        with self.assertRaisesRegex(InputError, "celestial WCS"):
            detect_positive_beam_scale_moment0_peaks(
                data,
                fits.Header({"BMAJ": 0.03, "BMIN": 0.02, "BPA": 0.0}),
                mask,
            )
        with self.assertRaisesRegex(InputError, "missing BMAJ, BMIN, BPA"):
            detect_positive_beam_scale_moment0_peaks(
                data, valid_header, mask
            )

    def test_rejected_label_is_removed_but_unrelated_label_remains(self):
        markers = np.zeros((8, 8), dtype=int)
        markers[1:3, 1:3] = 1
        markers[1:3, 4:6] = 2
        markers[5:7, 1:3] = 99
        audit = QTable(
            {
                "optical_label": [1, 2],
                "moment0_peak_supported": [True, False],
            }
        )

        filtered = remove_rejected_optical_regions(markers, audit)

        self.assertTrue(np.any(filtered == 1))
        self.assertFalse(np.any(filtered == 2))
        self.assertTrue(np.any(filtered == 99))
        self.assertTrue(np.any(markers == 2), "input array must not be mutated")

    def test_peak_filter_debug_products_record_exact_inputs_and_audit(self):
        _, header = make_wcs(width=6, height=6)
        header = add_beam(header)
        moment0 = np.arange(36, dtype=float).reshape(6, 6)
        peak_map = np.zeros((6, 6), dtype=np.uint8)
        peak_map[3, 4] = 1
        audit = QTable(
            {
                "Name": ["LS_DEBUG"],
                "type": ["SER"],
                "flux_g": [2.0],
                "optical_label": [5],
                "moment0_peak_supported": [True],
                "status": ["accepted"],
                "peak_x": [4.0],
                "peak_y": [3.0],
                "peak_ra": [10.0],
                "peak_dec": [-20.0],
                "peak_value": [22.0],
                "rejection_reason": [""],
            },
            units={"peak_ra": "deg", "peak_dec": "deg"},
        )

        with TemporaryDirectory() as temporary_directory:
            parent_name, peak_name, audit_name = (
                write_moment0_peak_filter_debug_products(
                    temporary_directory,
                    "36",
                    moment0,
                    header,
                    peak_map,
                    audit,
                )
            )
            written_audit = QTable.read(audit_name, format="ascii.ecsv")

            np.testing.assert_array_equal(fits.getdata(parent_name), moment0)
            np.testing.assert_array_equal(fits.getdata(peak_name), peak_map)
            self.assertEqual(written_audit["Name"][0], "LS_DEBUG")
            self.assertTrue(
                bool(written_audit["moment0_peak_supported"][0])
            )
            self.assertEqual(written_audit["status"][0], "accepted")

    def test_disabled_peak_filter_preserves_existing_dr10_marker_path(self):
        _, header = make_wcs()
        detected = np.ma.masked_array(
            np.array(
                [
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                ]
            ),
            mask=False,
        )
        selected = make_selected_table(
            WCS(header), [("LS_OLD_PATH", 1, 1, 1, 2.0, "SER")]
        )
        cfg = SimpleNamespace(
            input=SimpleNamespace(
                manual_input_tables=[None],
                manual_markers_only=False,
                auto_query_catalogue=True,
                galaxy_types=["SER"],
                filter_dr10_markers_by_moment0_peaks=False,
            ),
            general=SimpleNamespace(verbose=False, debug=False),
            internal=SimpleNamespace(
                auto_catalogue_path="automatic.csv", image_counter=0
            ),
            directories=SimpleNamespace(ancillary_directory="/unused"),
        )
        diagnostics = {}
        with (
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "detect_optical_sources",
                return_value=(detected, header, np.ones((8, 8))),
            ),
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "load_dr10_catalogue",
                return_value=QTable(),
            ),
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "select_dr10_counterparts",
                return_value=selected,
            ),
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "filter_dr10_counterparts_by_moment0_peaks"
            ) as peak_filter,
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "mask_source_from_table",
                return_value=detected,
            ),
        ):
            result, _ = set_optical_markers(
                cfg, "36", np.ones((3, 8, 8)), diagnostics=diagnostics
            )

        peak_filter.assert_not_called()
        self.assertIs(result, detected)
        self.assertIs(diagnostics["automatic_counterpart_table"], selected)

    def test_manual_catalogue_precedence_bypasses_peak_filter(self):
        _, header = make_wcs()
        detected = np.ma.masked_array(
            np.ones((8, 8), dtype=int), mask=False
        )
        manual_table = make_selected_table(
            WCS(header), [("MANUAL", 2, 2, 1, 1.0, "SER")]
        )
        cfg = SimpleNamespace(
            input=SimpleNamespace(
                manual_input_tables=["manual.csv"],
                manual_markers_only=False,
                auto_query_catalogue=True,
                galaxy_types=["SER"],
                filter_dr10_markers_by_moment0_peaks=True,
            ),
            general=SimpleNamespace(verbose=False, debug=False),
            internal=SimpleNamespace(
                auto_catalogue_path="automatic.csv", image_counter=0
            ),
            directories=SimpleNamespace(ancillary_directory="/unused"),
        )
        with (
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "detect_optical_sources",
                return_value=(detected, header, np.ones((8, 8))),
            ),
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "read_manual_table",
                return_value=manual_table,
            ),
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "load_dr10_catalogue"
            ) as automatic_loader,
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "filter_dr10_counterparts_by_moment0_peaks"
            ) as peak_filter,
            patch(
                "deblend_sofia_detections.deblending.deblending."
                "mask_source_from_table",
                return_value=detected,
            ),
        ):
            result, _ = set_optical_markers(
                cfg, "36", np.ones((3, 8, 8))
            )

        automatic_loader.assert_not_called()
        peak_filter.assert_not_called()
        self.assertIs(result, detected)

    def test_selected_dr10_row_is_used_as_positional_child_counterpart(self):
        with TemporaryDirectory() as temporary_directory:
            cubelet_directory = Path(temporary_directory) / "trial_cubelets"
            cubelet_directory.mkdir()
            header = fits.Header()
            header["BMAJ"] = 0.01
            header["CDELT1"] = -0.001
            fits.PrimaryHDU(
                np.zeros((2, 2, 2)), header=header
            ).writeto(cubelet_directory / "trial_1_cube.fits")
            source = QTable(
                {
                    "sofia_id": [1],
                    "sofia_ra": [10.0],
                    "sofia_dec": [-20.0],
                    "sofia_name": ["trial source"],
                },
                units={"sofia_ra": "deg", "sofia_dec": "deg"},
            )
            automatic_catalogue = QTable(
                {
                    "Name": ["LS_1000m200_42"],
                    "RA": [10.0001],
                    "DEC": [-20.0],
                    "type": ["SER"],
                    "flux_g": [10.0],
                },
                units={"RA": "deg", "DEC": "deg"},
            )
            cfg = SimpleNamespace(
                general=SimpleNamespace(
                    verbose=False, counterpart_region="Beam"
                )
            )

            result = search_counter_part(
                cfg,
                source,
                sofia_directory=temporary_directory,
                basename="trial",
                query="DR10",
                automatic_catalogue=automatic_catalogue,
            )

        self.assertTrue(bool(result["DR10_counterpart"][0]))
        self.assertEqual(result["DR10_Name"][0], "LS_1000m200_42")


if __name__ == "__main__":
    unittest.main()
