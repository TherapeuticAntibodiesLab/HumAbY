import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from humanize import (
    backmutation_counts,
    compare_backmutations,
    enrich_backmutation_report,
    format_germline_details,
    load_backmutation_proposals,
    parse_germline_metadata,
    write_candidate_fasta,
    write_candidate_summary,
)


class GermlineMetadataTests(unittest.TestCase):
    def test_loads_backmutation_proposals_from_fasta_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            artifact_dir = Path(temp_dir)
            (artifact_dir / "humanized_vh.fasta").write_text(
                ">vh1\nABCDE\n>vh2\nFGHIJ\n", encoding="utf-8"
            )
            (artifact_dir / "humanized_vl.fasta").write_text(
                ">vl1\nKLMNO\n", encoding="utf-8"
            )

            proposals = load_backmutation_proposals(artifact_dir)

        self.assertEqual(proposals, {"VH": ["ABCDE", "FGHIJ"], "VL": ["KLMNO"]})

    def test_counts_only_effective_backmutations_with_one_based_positions(self):
        changes = compare_backmutations("ABCDE", "ABXDY")

        self.assertEqual(changes, [
            {"position": 3, "from": "C", "to": "X"},
            {"position": 5, "from": "E", "to": "Y"},
        ])

    def test_rejects_backmutation_comparison_with_different_lengths(self):
        with self.assertRaisesRegex(ValueError, "different lengths"):
            compare_backmutations("ABC", "AB")

    def test_counts_backmutations_per_chain_and_candidate(self):
        raw = "Chain VH\nbackmt ABX\nbackmt DEF\nChain VL\nbackmt QRS\n"
        counts = backmutation_counts(raw, {"VH": ["ABC", "DEF"], "VL": ["QRT"]})

        self.assertEqual(counts, {"VH": [1, 0], "VL": [1], "total": 2})

    def test_parses_existing_database_title(self):
        title = (
            "VH_210771 Heavy chain framework | "
            "V:IGHV3-23*01 D:IGHD3-10*01 J:IGHJ4*02 | "
            "Human germline reconstruction"
        )

        metadata = parse_germline_metadata("ref|VH_210771|", title)

        self.assertEqual(metadata["chain"], "heavy")
        self.assertEqual(metadata["V"], "IGHV3-23*01")
        self.assertEqual(metadata["D"], "IGHD3-10*01")
        self.assertEqual(metadata["J"], "IGHJ4*02")

    def test_parses_enriched_database_title(self):
        title = (
            "VK_000231 Kappa light chain framework | "
            "V:IGKV1-39*01[accession=X59315,functionality=F] "
            "J:IGKJ1*01[accession=J00242,functionality=F] | "
            "Homo sapiens germline reconstruction with canonical FR4"
        )

        metadata = parse_germline_metadata("ref|VK_000231|", title)

        self.assertEqual(metadata["chain"], "kappa")
        self.assertEqual(metadata["V"], "IGKV1-39*01")
        self.assertEqual(metadata["V_accession"], "X59315")
        self.assertEqual(metadata["V_functionality"], "F")
        self.assertEqual(metadata["J"], "IGKJ1*01")
        self.assertEqual(metadata["J_accession"], "J00242")

    def test_writes_imgt_details_to_summary_and_fasta(self):
        candidate = {
            "source": "ref|VH_210771|",
            "identity": 62.4,
            "length": 4,
            "sequence": "QVQL",
            "germline": {
                "chain": "heavy",
                "species": "Homo sapiens",
                "V": "IGHV3-23*01",
                "D": "IGHD3-10*01",
                "J": "IGHJ4*02",
            },
        }

        summary = io.StringIO()
        write_candidate_summary([candidate], summary, "VH")
        self.assertIn("V=IGHV3-23*01", summary.getvalue())
        self.assertIn("D=IGHD3-10*01", summary.getvalue())
        self.assertIn("J=IGHJ4*02", summary.getvalue())

        with tempfile.TemporaryDirectory() as temp_dir:
            write_candidate_fasta([candidate], Path(temp_dir), "VH")
            fasta = (Path(temp_dir) / "humanized_vh.fasta").read_text()
        self.assertIn("V_IGHV3-23*01", fasta)
        self.assertIn("D_IGHD3-10*01", fasta)
        self.assertIn("J_IGHJ4*02", fasta)

        formatted = format_germline_details(candidate)
        self.assertIn("chain=heavy (IGH)", formatted)
        self.assertIn("species=Homo sapiens", formatted)

    @patch("humanize.find_homologous_frameworks")
    def test_enriches_backmutation_report_from_backmutated_sequences(self, find_frameworks):
        def fake_find(sequence, chain_type, database_path):
            self.assertEqual(database_path, "human_db")
            if sequence == "BACKMUTATEDVH":
                return [
                    {
                        "seq_id": "ref|VH_BACK|",
                        "identity": 77.7,
                        "sequence": sequence,
                        "evalue": 1e-20,
                        "bitscore": 120.0,
                        "germline": {"chain": "heavy", "V": "IGHV7-7*01", "D": "IGHD7-7*01", "J": "IGHJ7*01"},
                    },
                    {
                        "seq_id": "ref|VH_SECOND|",
                        "identity": 75.5,
                        "sequence": sequence,
                        "evalue": 2e-18,
                        "bitscore": 110.0,
                        "germline": {"chain": "heavy", "V": "IGHV6-6*01", "D": "IGHD6-6*01", "J": "IGHJ6*01"},
                    },
                ]
            if sequence == "BACKMUTATEDVL":
                return [{
                    "seq_id": "ref|VL_BACK|",
                    "identity": 88.8,
                    "sequence": sequence,
                    "evalue": 1e-30,
                    "bitscore": 130.0,
                    "germline": {"chain": "kappa", "V": "IGKV8-8*01", "J": "IGKJ8*01"},
                }]
            return []

        find_frameworks.side_effect = fake_find
        raw_report = (
            "\n Chain VH\n"
            "hGerm copied-from-opt1 germline\n"
            "backmt BACKMUTATEDVH\n"
            "\n Chain VL\n"
            "hGerm copied-from-opt1 light germline\n"
            "backmt BACKMUTATEDVL\n"
        )

        backmutation_log = io.StringIO()
        enriched = enrich_backmutation_report(
            raw_report,
            "human_db",
            backmutation_log,
            {"VH": ["BACKMUTATEDVA"], "VL": ["BACKMUTATEDVL"]},
        )
        log_text = backmutation_log.getvalue()

        self.assertIn("BACKMUTATION RESULTS", enriched)
        self.assertIn("CHAIN VH", enriched)
        self.assertIn("VH candidate 1", enriched)
        self.assertIn("Human germline template used for grafting:", enriched)
        self.assertIn("Backmutated sequence:", enriched)
        self.assertIn("Back-mutations realizadas: 1", enriched)
        self.assertIn("13: A -> H", enriched)
        self.assertIn("Germline determined from this backmt sequence:", enriched)
        self.assertIn("Source: ref|VH_BACK|", enriched)
        self.assertIn("Identity: 77.7%", enriched)
        self.assertIn("V=IGHV7-7*01", enriched)
        self.assertIn("D=IGHD7-7*01", enriched)
        self.assertIn("J=IGHJ7*01", enriched)
        self.assertIn("CHAIN VL", enriched)
        self.assertIn("VL candidate 1", enriched)
        self.assertIn("Source: ref|VL_BACK|", enriched)
        self.assertIn("V=IGKV8-8*01", enriched)
        self.assertIn("J=IGKJ8*01", enriched)
        self.assertNotIn("VH_candidate", enriched)
        self.assertLess(enriched.index("backmt: BACKMUTATEDVH"), enriched.index("Source: ref|VH_BACK|"))
        self.assertIn("Found VH backmt sequence #1", log_text)
        self.assertIn("Analyzing VH backmutated sequence against human germline database", log_text)
        self.assertIn("Germline decision method: run blastp", log_text)
        self.assertIn("Human BLAST search parameters:", log_text)
        self.assertIn("Human germline candidates returned for VH: 2; ranking criterion=bitscore descending", log_text)
        self.assertIn("Candidate rank 1: source=ref|VH_BACK|", log_text)
        self.assertIn("Candidate rank 2: source=ref|VH_SECOND|", log_text)
        self.assertIn("Best VH backmutated germline match: source=ref|VH_BACK|", log_text)
        self.assertIn("Selection reason: rank #1 by highest BLAST bitscore (120.0)", log_text)


if __name__ == "__main__":
    unittest.main()
