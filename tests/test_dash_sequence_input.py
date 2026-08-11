import unittest
from unittest.mock import patch

import dash_app
from dash_app import (
    launch_run,
    normalize_sequence_text,
    parse_backmutation_display,
    parse_candidate_summary,
    sequence_pair_to_input_text,
    validate_sequence_pair,
)


class DashSequenceInputTests(unittest.TestCase):
    def test_parses_candidates_for_visual_results(self):
        summary = (
            "VH Chain: 1 candidates\n"
            "  1. ref|VH_001| (identity: 72.5%, length: 117 AA)\n"
            "     IMGT germline: chain=heavy; V=IGHV1-2*01; J=IGHJ4*01\n"
        )

        self.assertEqual(parse_candidate_summary(summary), [{
            "chain": "VH", "candidate": 1, "source": "ref|VH_001|",
            "identity": 72.5, "length": 117,
            "germline": "chain=heavy; V=IGHV1-2*01; J=IGHJ4*01",
        }])

    def test_parses_backmutations_for_visual_results(self):
        report = (
            "VH candidate 1\n"
            "Back-mutations performed: 2\n"
            "Changes (position: humanized -> backmutated):\n"
            "  48: L -> I\n"
            "  71: R -> A\n"
        )

        self.assertEqual(parse_backmutation_display(report), {
            ("VH", 1): {
                "count": 2,
                "changes": [
                    {"position": 48, "from": "L", "to": "I"},
                    {"position": 71, "from": "R", "to": "A"},
                ],
            }
        })

    def test_normalizes_wrapped_sequence(self):
        self.assertEqual(normalize_sequence_text(" qvql\n vqsg "), "QVQLVQSG")

    def test_validates_pasted_pair(self):
        vh, vl = validate_sequence_pair("A CDEFGHIKLMNPQRSTVWY" * 3, "V\nACDEFGHIKLMNPQRSTWY" * 3)
        self.assertGreaterEqual(len(vh), 50)
        self.assertGreaterEqual(len(vl), 50)
        self.assertNotIn(" ", vh)
        self.assertNotIn("\n", vl)

    def test_rejects_invalid_amino_acid(self):
        with self.assertRaisesRegex(ValueError, "VH contains invalid"):
            validate_sequence_pair("B" * 60, "A" * 60)

    def test_creates_pipeline_input(self):
        self.assertEqual(sequence_pair_to_input_text("VHSEQ", "VLSEQ"), "VHSEQ\nVLSEQ\n")

    @patch.object(dash_app, "start_background_run")
    @patch.object(dash_app, "create_run_structure")
    def test_launch_uses_pasted_sequences_before_upload(self, create_run, start_run):
        create_run.return_value = {"run_id": "test-run"}
        vh = "ACDEFGHIKLMNPQRSTVWY" * 3
        vl = "YWVTSRQPNMLKIHGFEDCA" * 3

        _, store = launch_run(
            1,
            vh,
            vl,
            "ignored-upload",
            "ignored.txt",
            "pasted-input",
            "graft",
            None,
        )

        create_run.assert_called_once_with(
            "pasted-input",
            f"{vh}\n{vl}\n",
            "graft",
            None,
            False,
        )
        start_run.assert_called_once_with("test-run")
        self.assertEqual(store, {"run_id": "test-run"})

    @patch.object(dash_app, "start_background_run")
    @patch.object(dash_app, "create_run_structure")
    def test_launch_passes_graft_backmutation_mode(self, create_run, start_run):
        create_run.return_value = {"run_id": "backmutation-run"}
        vh = "ACDEFGHIKLMNPQRSTVWY" * 3
        vl = "YWVTSRQPNMLKIHGFEDCA" * 3

        _, store = launch_run(
            1,
            vh,
            vl,
            None,
            None,
            "backmutation-input",
            "graft_backmutation",
            None,
        )

        create_run.assert_called_once_with(
            "backmutation-input",
            f"{vh}\n{vl}\n",
            "graft_backmutation",
            None,
            True,
        )
        start_run.assert_called_once_with("backmutation-run")
        self.assertEqual(store, {"run_id": "backmutation-run"})


if __name__ == "__main__":
    unittest.main()
