import unittest
from contextlib import redirect_stdout
from io import StringIO

import numpy as np

from deblend_sofia_detections.deblending.debug_visualization import (
    catalogue_positions_from_table,
    contour_levels,
    deduplicate_positions,
    marker_centroids,
    matched_counterpart_positions,
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
