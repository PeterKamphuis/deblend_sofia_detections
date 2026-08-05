import os
import tempfile
import unittest
from contextlib import chdir

from omegaconf import OmegaConf

from deblend_sofia_detections.config.config import defaults
from deblend_sofia_detections.config.functions import setup_config


class PrintedConfigurationTests(unittest.TestCase):
    def test_printed_example_uses_recommended_deblending_settings(self):
        with tempfile.TemporaryDirectory() as directory, chdir(directory):
            with self.assertRaises(SystemExit):
                setup_config(["print_examples=true"])

            generated = OmegaConf.load(
                os.path.join(
                    directory, "deblend_sofia_detections_default.yml"
                )
            )
            self.assertTrue(generated.input.auto_query_catalogue)
            self.assertTrue(
                generated.input.filter_dr10_markers_by_moment0_peaks
            )
            self.assertTrue(
                generated.input[
                    "deblend_optical_regions_with_multiple_moment0_peaks"
                ]
            )
            self.assertFalse(generated.input.use_peak_deblending)
            self.assertTrue(os.path.isfile("sofia_template.par"))

    def test_runtime_compatibility_defaults_are_unchanged(self):
        runtime = OmegaConf.structured(defaults)

        self.assertFalse(runtime.input.auto_query_catalogue)
        self.assertFalse(
            runtime.input.filter_dr10_markers_by_moment0_peaks
        )
        self.assertFalse(
            runtime.input.deblend_optical_regions_with_multiple_moment0_peaks
        )
        self.assertTrue(runtime.input.use_peak_deblending)


if __name__ == "__main__":
    unittest.main()
