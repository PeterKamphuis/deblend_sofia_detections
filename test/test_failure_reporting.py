import tempfile
import unittest
from pathlib import Path

from deblend_sofia_detections.support.failure_reporting import (
    SourceFailure,
    format_failure_report,
    write_failure_report,
)


class FailureReportingTests(unittest.TestCase):
    def setUp(self):
        self.failure = SourceFailure(
            source_id="1",
            cube_name="/data/field_1_cube.fits",
            exception_type="AttributeError",
            reason="'NoneType' object has no attribute 'wcs'",
            traceback_text="Traceback (most recent call last):\nexample",
        )

    def test_report_contains_failure_reason_and_traceback(self):
        report = format_failure_report(
            [self.failure],
            requested_count=9,
            succeeded_count=8,
            generated_at="2026-07-23T12:00:00+02:00",
        )

        self.assertIn("Requested sources: 9", report)
        self.assertIn("Succeeded: 8", report)
        self.assertIn("Failed: 1", report)
        self.assertIn("Source ID: 1", report)
        self.assertIn("Exception: AttributeError", report)
        self.assertIn(self.failure.reason, report)
        self.assertIn(self.failure.traceback_text, report)

    def test_successful_run_replaces_stale_failures(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            report_path = Path(temporary_directory) / "deblend_failures.log"
            write_failure_report(
                report_path,
                [self.failure],
                requested_count=1,
                succeeded_count=0,
            )
            write_failure_report(
                report_path,
                [],
                requested_count=2,
                succeeded_count=2,
            )

            report = report_path.read_text(encoding="utf-8")

        self.assertIn("Failed: 0", report)
        self.assertIn("No source-level failures were recorded.", report)
        self.assertNotIn("Source ID: 1", report)


if __name__ == "__main__":
    unittest.main()
