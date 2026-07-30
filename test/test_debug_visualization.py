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
    trial_hi_components_from_table,
    write_hi_component_debug_overlay,
    write_hi_component_debug_overlay_safely,
    write_source_debug_overlay_safely,
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
