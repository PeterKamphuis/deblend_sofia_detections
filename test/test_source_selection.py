import unittest

from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.support.source_selection import select_source_ids


class SelectSourceIdsTests(unittest.TestCase):
    def test_empty_allowlist_selects_every_catalogue_source(self):
        catalogue_ids = [1, 2, 3]

        self.assertEqual(select_source_ids(catalogue_ids, []), catalogue_ids)

    def test_selection_preserves_catalogue_types_and_order(self):
        catalogue_ids = [3, 1, 2]

        selected = select_source_ids(catalogue_ids, ["2", "3"])

        self.assertEqual(selected, [3, 2])
        self.assertTrue(all(isinstance(source_id, int) for source_id in selected))

    def test_duplicate_requested_ids_are_processed_once(self):
        selected = select_source_ids([1, 2, 3], ["2", "2"])

        self.assertEqual(selected, [2])

    def test_missing_source_ids_raise_an_input_error(self):
        with self.assertRaisesRegex(InputError, "4.*9"):
            select_source_ids([1, 2, 3], ["2", "9", "4"])


if __name__ == "__main__":
    unittest.main()
