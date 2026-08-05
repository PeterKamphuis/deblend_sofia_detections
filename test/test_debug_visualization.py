import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from deblend_sofia_detections.deblending.debug_visualization import (
    catalogue_positions_from_table,
    component_colour_mapping,
    contour_levels,
    deduplicate_positions,
    marker_centroids,
    matched_counterpart_positions,
    pv_projection_data,
    trial_hi_components_from_table,
    write_hi_component_debug_overlay,
    write_hi_component_debug_overlay_safely,
    write_hi_component_pv_debug_overlays,
    write_source_debug_overlay_safely,
    write_source_pv_debug_overlays,
)


class FakeTable(list):
    def __init__(self, rows):
        super().__init__(rows)
        self.colnames = list(rows[0]) if rows else []

    def __getitem__(self, key):
        if isinstance(key, str):
            return np.asarray([row[key] for row in self])
        return super().__getitem__(key)


class DebugVisualizationTests(unittest.TestCase):
    @staticmethod
    def make_cube_wcs(width=20, height=16, channels=12):
        from astropy.wcs import WCS

        wcs = WCS(naxis=3)
        wcs.wcs.crpix = [width / 2.0, height / 2.0, 1.0]
        wcs.wcs.cdelt = np.array([-0.001, 0.001, 20000.0])
        wcs.wcs.crval = [10.0, -20.0, 10000000.0]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "VRAD"]
        wcs.wcs.cunit = ["deg", "deg", "m/s"]
        header = wcs.to_header()
        return wcs, header

    def test_marker_centroids_are_reported_for_each_positive_label(self):
        markers = np.array(
            [
                [0, 1, 1, 0],
                [0, 1, 1, 0],
                [2, 0, 0, 0],
            ]
        )

        self.assertEqual(
            marker_centroids(markers),
            [(1.5, 0.5, 1), (0.0, 2.0, 2)],
        )

    def test_contours_only_use_values_inside_the_source_mask(self):
        data = np.array([[1.0, 2.0], [1000.0, 4.0]])
        mask = np.array([[True, True], [False, True]])

        levels = contour_levels(data, mask)

        self.assertEqual(len(levels), 3)
        self.assertTrue(all(1.0 < level < 4.0 for level in levels))
        self.assertEqual(levels, sorted(levels))

    def test_pv_projection_shapes_and_velocity_units(self):
        _, header = self.make_cube_wcs(width=10, height=8, channels=6)
        cube = np.arange(6 * 8 * 10, dtype=float).reshape(6, 8, 10)
        mask = np.zeros_like(cube, dtype=int)
        mask[1:5, 2:7, 3:9] = 1

        ra_projection = pv_projection_data(cube, header, mask, "ra")
        dec_projection = pv_projection_data(cube, header, mask, "dec")

        self.assertEqual(ra_projection["background"].shape, (6, 10))
        self.assertEqual(ra_projection["source"].shape, (6, 10))
        self.assertEqual(dec_projection["background"].shape, (6, 8))
        self.assertEqual(dec_projection["source"].shape, (6, 8))
        self.assertAlmostEqual(ra_projection["velocity_kms"][0], 10000.0)
        self.assertAlmostEqual(ra_projection["velocity_kms"][1], 10020.0)
        self.assertFalse(bool(ra_projection["support"][0, 4]))
        self.assertTrue(bool(ra_projection["support"][2, 4]))

    def test_frequency_pv_axis_uses_rest_frequency_for_radio_velocity(self):
        from astropy.wcs import WCS

        rest_frequency = 1.420405751e9
        wcs = WCS(naxis=3)
        wcs.wcs.crpix = [3.0, 3.0, 1.0]
        wcs.wcs.cdelt = np.array([-0.001, 0.001, -100000.0])
        wcs.wcs.crval = [10.0, -20.0, rest_frequency]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "FREQ"]
        wcs.wcs.cunit = ["deg", "deg", "Hz"]
        header = wcs.to_header()
        header["RESTFRQ"] = rest_frequency
        cube = np.ones((4, 6, 6), dtype=float)

        projection = pv_projection_data(cube, header, cube > 0, "ra")

        self.assertAlmostEqual(projection["velocity_kms"][0], 0.0)
        self.assertTrue(np.all(np.diff(projection["velocity_kms"]) > 0))

    def test_component_colours_follow_sorted_ids_and_are_distinct(self):
        colours = component_colour_mapping(["10", 2, "1"])

        self.assertEqual(list(colours), ["1", "2", "10"])
        self.assertEqual(len(set(colours.values())), 3)

    def test_trial_components_skip_missing_maps_and_keep_invalid_centres(self):
        from astropy.io import fits
        from astropy.table import Table

        with TemporaryDirectory() as temporary_directory:
            cubelet_directory = Path(temporary_directory)
            fits.PrimaryHDU(np.ones((4, 4))).writeto(
                cubelet_directory / "trial_1_mom0.fits"
            )
            fits.PrimaryHDU(np.ones((4, 4))).writeto(
                cubelet_directory / "trial_2_mom0.fits"
            )
            table = Table(
                {
                    "id": [2, 3, 1],
                    "ra": [10.2, 10.3, np.nan],
                    "dec": [-20.2, -20.3, np.nan],
                }
            )
            output = StringIO()
            with redirect_stdout(output):
                components = trial_hi_components_from_table(
                    table,
                    cubelet_directory=cubelet_directory,
                    basename="trial",
                )

        self.assertEqual(
            [component["id"] for component in components],
            ["1", "2"],
        )
        self.assertTrue(np.isnan(components[0]["ra_deg"]))
        self.assertEqual(components[1]["ra_deg"], 10.2)
        self.assertIn("component 3", output.getvalue())

    def test_two_trial_children_render_with_catalogue_centres(self):
        from astropy.io import fits
        from astropy.wcs import WCS

        with TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            wcs = WCS(naxis=2)
            wcs.wcs.crpix = [20.0, 20.0]
            wcs.wcs.cdelt = np.array([-0.001, 0.001])
            wcs.wcs.crval = [10.0, -20.0]
            wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            header = wcs.to_header()

            optical_name = directory / "optical.fits"
            fits.PrimaryHDU(
                np.arange(1600, dtype=float).reshape(40, 40),
                header=header,
            ).writeto(optical_name)

            child_names = []
            centres = []
            yy, xx = np.indices((40, 40))
            for component_id, x_centre, y_centre in (
                ("1", 14, 18),
                ("2", 26, 22),
            ):
                child_data = np.exp(
                    -(
                        (xx - x_centre) ** 2
                        + (yy - y_centre) ** 2
                    )
                    / 12.0
                )
                child_name = directory / f"trial_{component_id}_mom0.fits"
                fits.PrimaryHDU(child_data, header=header).writeto(child_name)
                child_names.append(child_name)
                centres.append(
                    wcs.pixel_to_world_values(x_centre, y_centre)
                )

            output_name = directory / "component_overlay.png"
            output_path = write_hi_component_debug_overlay(
                optical_image_name=optical_name,
                moment0_data=np.arange(
                    1600, dtype=float
                ).reshape(40, 40),
                moment0_header=header,
                source_mask=np.ones((3, 40, 40)),
                marker_data=np.zeros((40, 40)),
                catalogue_positions=[],
                components=[
                    {
                        "id": str(index + 1),
                        "moment0_name": str(child_name),
                        "ra_deg": centres[index][0],
                        "dec_deg": centres[index][1],
                    }
                    for index, child_name in enumerate(child_names)
                ],
                output_name=output_name,
                source_id="36",
            )

            self.assertEqual(output_path, output_name)
            self.assertTrue(output_name.is_file())
            self.assertGreater(output_name.stat().st_size, 0)

    def test_component_safe_writer_ignores_empty_component_list(self):
        result = write_hi_component_debug_overlay_safely(
            components=[],
            source_id="36",
        )

        self.assertIsNone(result)

    def test_catalogue_and_component_pv_overlays_render_both_axes(self):
        from astropy.io import fits

        with TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            cube_wcs, cube_header = self.make_cube_wcs()
            channels, height, width = 12, 16, 20
            zz, yy, xx = np.indices((channels, height, width))
            cube = np.exp(
                -(
                    (xx - 9.0) ** 2 / 15.0
                    + (yy - 8.0) ** 2 / 10.0
                    + (zz - 5.0) ** 2 / 5.0
                )
            )
            mask = (cube > 0.08).astype(np.int16)
            optical_header = cube_wcs.celestial.to_header()
            optical_name = directory / "optical.fits"
            fits.PrimaryHDU(
                np.sum(cube, axis=0), header=optical_header
            ).writeto(optical_name)
            markers = np.zeros((height, width), dtype=int)
            markers[4:12, 5:14] = 1
            ra_deg, dec_deg = cube_wcs.celestial.pixel_to_world_values(9, 8)
            positions = [
                {
                    "ra_deg": ra_deg,
                    "dec_deg": dec_deg,
                    "name": "LS_TEST",
                    "catalogue": "Legacy Surveys DR10",
                    "marker_status": "accepted",
                }
            ]

            catalogue_ra = directory / "catalogue_ra.png"
            catalogue_dec = directory / "catalogue_dec.png"
            catalogue_paths = write_source_pv_debug_overlays(
                cube_data=cube,
                cube_header=cube_header,
                source_mask=mask,
                optical_image_name=optical_name,
                marker_data=markers,
                catalogue_positions=positions,
                output_ra_name=catalogue_ra,
                output_dec_name=catalogue_dec,
                source_id="36",
            )

            child_cube_name = directory / "child_cube.fits"
            child_mask_name = directory / "child_mask.fits"
            fits.PrimaryHDU(cube, header=cube_header).writeto(child_cube_name)
            fits.PrimaryHDU(mask, header=cube_header).writeto(child_mask_name)
            component_ra = directory / "component_ra.png"
            component_dec = directory / "component_dec.png"
            component_paths = write_hi_component_pv_debug_overlays(
                cube_data=cube,
                cube_header=cube_header,
                source_mask=mask,
                optical_image_name=optical_name,
                marker_data=markers,
                catalogue_positions=positions,
                components=[
                    {
                        "id": "1",
                        "cube_name": str(child_cube_name),
                        "mask_name": str(child_mask_name),
                        "ra_deg": ra_deg,
                        "dec_deg": dec_deg,
                        "velocity_kms": 10100.0,
                    }
                ],
                output_ra_name=component_ra,
                output_dec_name=component_dec,
                source_id="36",
            )

            self.assertEqual(catalogue_paths, (catalogue_ra, catalogue_dec))
            self.assertEqual(component_paths, (component_ra, component_dec))
            for output_name in (
                catalogue_ra, catalogue_dec, component_ra, component_dec
            ):
                self.assertTrue(output_name.is_file())
                self.assertGreater(output_name.stat().st_size, 0)

    def test_source_overlay_renders_accepted_and_rejected_dr10_positions(self):
        from astropy.io import fits
        from astropy.table import Table
        from astropy.wcs import WCS

        with TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            wcs = WCS(naxis=2)
            wcs.wcs.crpix = [10.0, 10.0]
            wcs.wcs.cdelt = np.array([-0.001, 0.001])
            wcs.wcs.crval = [10.0, -20.0]
            wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            header = wcs.to_header()
            optical_name = directory / "optical.fits"
            fits.PrimaryHDU(
                np.arange(400, dtype=float).reshape(20, 20),
                header=header,
            ).writeto(optical_name)
            accepted_ra, accepted_dec = wcs.pixel_to_world_values(8, 10)
            rejected_ra, rejected_dec = wcs.pixel_to_world_values(12, 10)
            output_name = directory / "overlay.png"
            catalogue_name = directory / "positions.ecsv"

            result = write_source_debug_overlay_safely(
                optical_image_name=optical_name,
                moment0_data=np.arange(400, dtype=float).reshape(20, 20),
                moment0_header=header,
                source_mask=np.ones((3, 20, 20)),
                marker_data=np.zeros((20, 20)),
                catalogue_positions=[
                    {
                        "ra_deg": accepted_ra,
                        "dec_deg": accepted_dec,
                        "name": "LS_ACCEPTED",
                        "catalogue": "Legacy Surveys DR10",
                        "marker_status": "accepted",
                    },
                    {
                        "ra_deg": rejected_ra,
                        "dec_deg": rejected_dec,
                        "name": "LS_REJECTED",
                        "catalogue": "Legacy Surveys DR10",
                        "marker_status": "rejected",
                    },
                ],
                output_name=output_name,
                catalogue_output_name=catalogue_name,
                source_id="36",
                marker_mode="automatic + DR10 + moment-0 peak filter",
            )

            self.assertEqual(result, output_name)
            self.assertTrue(output_name.is_file())
            self.assertGreater(output_name.stat().st_size, 0)
            self.assertEqual(len(Table.read(catalogue_name)), 2)

    def test_invalid_child_wcs_does_not_hide_valid_child(self):
        from astropy.io import fits
        from astropy.wcs import WCS

        with TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            wcs = WCS(naxis=2)
            wcs.wcs.crpix = [10.0, 10.0]
            wcs.wcs.cdelt = np.array([-0.001, 0.001])
            wcs.wcs.crval = [10.0, -20.0]
            wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            header = wcs.to_header()
            yy, xx = np.indices((20, 20))
            data = np.exp(-((xx - 10) ** 2 + (yy - 10) ** 2) / 8.0)

            optical_name = directory / "optical.fits"
            valid_name = directory / "trial_1_mom0.fits"
            invalid_name = directory / "trial_2_mom0.fits"
            fits.PrimaryHDU(data, header=header).writeto(optical_name)
            fits.PrimaryHDU(data, header=header).writeto(valid_name)
            fits.PrimaryHDU(data).writeto(invalid_name)

            output_name = directory / "component_overlay.png"
            output = StringIO()
            with redirect_stdout(output):
                result = write_hi_component_debug_overlay_safely(
                    optical_image_name=optical_name,
                    moment0_data=data,
                    moment0_header=header,
                    source_mask=np.ones((2, 20, 20)),
                    marker_data=np.zeros((20, 20)),
                    catalogue_positions=[],
                    components=[
                        {
                            "id": "1",
                            "moment0_name": str(valid_name),
                            "ra_deg": 10.0,
                            "dec_deg": -20.0,
                        },
                        {
                            "id": "2",
                            "moment0_name": str(invalid_name),
                            "ra_deg": np.nan,
                            "dec_deg": np.nan,
                        },
                    ],
                    output_name=output_name,
                    source_id="36",
                )

            self.assertEqual(result, output_name)
            self.assertTrue(output_name.is_file())
            self.assertIn("component 2", output.getvalue())

    def test_later_catalogue_match_replaces_duplicate_input_position(self):
        positions = [
            {
                "ra_deg": 10.0,
                "dec_deg": -20.0,
                "name": "Input A",
                "catalogue": "Manual input",
            },
            {
                "ra_deg": 10.0 + 1e-8,
                "dec_deg": -20.0,
                "name": "Matched A",
                "catalogue": "Manual match",
            },
        ]

        unique = deduplicate_positions(positions)

        self.assertEqual(len(unique), 1)
        self.assertEqual(unique[0]["name"], "Matched A")
        self.assertEqual(unique[0]["catalogue"], "Manual match")

    def test_manual_catalogue_positions_are_extracted_case_insensitively(self):
        table = FakeTable(
            [{"name": "Galaxy A", "ra": 10.5, "dec": -20.25}]
        )

        positions = catalogue_positions_from_table(table)

        self.assertEqual(
            positions,
            [
                {
                    "ra_deg": 10.5,
                    "dec_deg": -20.25,
                    "name": "Galaxy A",
                    "catalogue": "Manual input",
                }
            ],
        )

    def test_peak_filter_audit_positions_keep_accepted_rejected_status(self):
        table = FakeTable(
            [
                {
                    "Name": "LS_ACCEPTED",
                    "RA": 10.5,
                    "DEC": -20.25,
                    "moment0_peak_supported": True,
                },
                {
                    "Name": "LS_REJECTED",
                    "RA": 10.6,
                    "DEC": -20.35,
                    "moment0_peak_supported": False,
                },
            ]
        )

        positions = catalogue_positions_from_table(
            table, catalogue="Legacy Surveys DR10"
        )

        self.assertEqual(
            [position["marker_status"] for position in positions],
            ["accepted", "rejected"],
        )
        self.assertEqual(
            deduplicate_positions(positions)[1]["marker_status"],
            "rejected",
        )

    def test_confirmed_manual_match_is_preferred_over_ned_match(self):
        table = FakeTable(
            [
                {
                    "Manual_spectroscopic": True,
                    "Manual_RA": 10.0,
                    "Manual_DEC": -20.0,
                    "Manual_Name": "Manual A",
                    "NED_spectroscopic": True,
                    "NED_RA": 11.0,
                    "NED_DEC": -21.0,
                    "NED_Object Name": "NED A",
                }
            ]
        )

        positions = matched_counterpart_positions(table)

        self.assertEqual(len(positions), 1)
        self.assertEqual(positions[0]["name"], "Manual A")
        self.assertEqual(positions[0]["catalogue"], "Manual match")

    def test_dr10_match_is_reported_when_no_manual_match_exists(self):
        table = FakeTable(
            [
                {
                    "Manual_spectroscopic": False,
                    "DR10_counterpart": True,
                    "DR10_RA": 10.25,
                    "DR10_DEC": -20.5,
                    "DR10_Name": "LS_1000m200_42",
                    "NED_spectroscopic": False,
                }
            ]
        )

        positions = matched_counterpart_positions(table)

        self.assertEqual(len(positions), 1)
        self.assertEqual(positions[0]["name"], "LS_1000m200_42")
        self.assertEqual(
            positions[0]["catalogue"], "Legacy Surveys DR10 match"
        )

    def test_plotting_failure_is_reported_without_raising(self):
        output = StringIO()
        with redirect_stdout(output):
            result = write_source_debug_overlay_safely(
                optical_image_name="/missing/optical.fits",
                moment0_data=np.zeros((2, 2)),
                moment0_header={},
                source_mask=np.ones((2, 2)),
                marker_data=np.zeros((2, 2)),
                catalogue_positions=[],
                output_name="/missing/overlay.png",
                catalogue_output_name="/missing/positions.ecsv",
                source_id="36",
            )

        self.assertIsNone(result)
        self.assertIn("SoFiA source 36", output.getvalue())


if __name__ == "__main__":
    unittest.main()
