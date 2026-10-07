"""Synthetic, offline invariants for the campaign's submission helper."""

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "skills/tooluniverse-protein-design-campaign/scripts/preflight_submission.py"
)
spec = importlib.util.spec_from_file_location("campaign_preflight", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PreflightTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.csv = self.root / "submission.csv"
        self.args = SimpleNamespace(
            csv=self.csv,
            methods=None,
            source_csv=None,
            min_length=10,
            max_length=250,
            max_designs=20,
            max_bytes=5_000_000,
            classes=module.DEFAULT_CLASSES,
            paired_classes=module.DEFAULT_PAIRED_CLASSES,
            default_class="single_chain",
        )

    def write_rows(self, rows, path=None):
        with (path or self.csv).open("w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["name", "sequence", "molecule_class"]
            )
            writer.writeheader()
            writer.writerows(rows)

    def row(self, name="synthetic_a", sequence="ACDEFGHIKL", cls="single_chain"):
        return {"name": name, "sequence": sequence, "molecule_class": cls}

    def check(self):
        return module.preflight(self.args)

    def test_valid_hash_and_no_input_mutation(self):
        self.write_rows([self.row()])
        before = self.csv.read_bytes()
        result = self.check()
        self.assertTrue(result["format_and_disclosure_checks_pass"])
        self.assertEqual(result["candidate_count"], 1)
        self.assertEqual(result["candidates"][0]["chain_lengths"], [10])
        self.assertEqual(self.csv.read_bytes(), before)
        self.assertIn("official_novelty", result["not_checked"])

    def test_paired_lengths_are_per_chain(self):
        self.write_rows(
            [self.row(sequence="A" * 250 + ":" + "G" * 250, cls="fab_kappa")]
        )
        self.assertTrue(self.check()["format_and_disclosure_checks_pass"])
        self.write_rows(
            [self.row(sequence="A" * 251 + ":" + "G" * 250, cls="fab_kappa")]
        )
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])

    def test_duplicate_name_and_sequence(self):
        self.write_rows([self.row(), self.row()])
        issues = {i["category"] for i in self.check()["issues"]}
        self.assertTrue({"duplicate_name", "duplicate_sequence"}.issubset(issues))

    def test_source_mismatch_including_class(self):
        self.write_rows([self.row()])
        source = self.root / "source.csv"
        self.write_rows([self.row(cls="nanobody")], source)
        self.args.source_csv = source
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.write_rows([self.row()], source)
        self.assertTrue(self.check()["format_and_disclosure_checks_pass"])

    def test_missing_source_entry(self):
        self.write_rows([self.row()])
        source = self.root / "source.csv"
        self.write_rows([self.row(name="synthetic_b")], source)
        self.args.source_csv = source
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])

    def test_formula_and_control_characters_in_name(self):
        self.write_rows([self.row(name="  =1+1")])
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.write_rows([self.row(name="synthetic\nname")])
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])

    def test_secret_is_flagged_without_echo(self):
        self.write_rows([self.row()])
        methods = self.root / "methods.txt"
        secret = "SYNTHETIC_SECRET_NEVER_REAL"
        methods.write_text("https://example.org/result?token=" + secret)
        self.args.methods = methods
        result = self.check()
        self.assertFalse(result["format_and_disclosure_checks_pass"])
        self.assertNotIn(secret, json.dumps(result))
        self.assertNotIn("synthetic_a", json.dumps(result))
        self.assertNotIn("ACDEFGHIKL", json.dumps(result))

    def test_paths_email_and_private_ip(self):
        self.write_rows([self.row()])
        methods = self.root / "methods.txt"
        methods.write_text("/home/example/run example@example.org 172.16.1.2")
        self.args.methods = methods
        issues = {i["category"] for i in self.check()["issues"]}
        self.assertTrue(
            {"private_path", "email_address", "private_ip_or_localhost"}.issubset(
                issues
            )
        )

    def test_portal_vocabulary_can_be_configured(self):
        self.write_rows([self.row(cls="protein")])
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.args.classes = "protein,nanobody"
        self.assertTrue(self.check()["format_and_disclosure_checks_pass"])

    def test_scfv_chain_convention_can_be_configured(self):
        self.write_rows([self.row(cls="scfv")])
        self.assertTrue(self.check()["format_and_disclosure_checks_pass"])
        self.args.paired_classes = "scfv,fab_kappa,fab_lambda"
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.write_rows([self.row(sequence="ACDEFGHIKL:GGGGGGGGGG", cls="scfv")])
        self.assertTrue(self.check()["format_and_disclosure_checks_pass"])

    def test_custom_paired_class(self):
        self.write_rows([self.row(sequence="ACDEFGHIKL:GGGGGGGGGG", cls="heterodimer")])
        self.args.classes = "heterodimer"
        self.args.paired_classes = "heterodimer"
        self.assertTrue(self.check()["format_and_disclosure_checks_pass"])

    def test_count_invalid_sequence_and_wrong_chain_format(self):
        self.write_rows(
            [self.row(), self.row(name="synthetic_b", sequence="GGGGGGGGGG")]
        )
        self.args.max_designs = 1
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.args.max_designs = 20
        self.write_rows([self.row(sequence="ACDEFGHIKX")])
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.write_rows([self.row(sequence="ACDEFGHIKL:GGGGGGGGGG")])
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])

    def test_duplicate_headers_and_extra_columns(self):
        self.csv.write_text("name,sequence,name\na,ACDEFGHIKL,b\n")
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.csv.write_text("name,sequence,internal_log\na,ACDEFGHIKL,unreviewed\n")
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])

    def test_bad_bytes_and_limits(self):
        self.csv.write_bytes(b"\xff\xfe")
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.write_rows([self.row()])
        self.args.max_bytes = 1
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])
        self.args.min_length = 251
        self.assertFalse(self.check()["format_and_disclosure_checks_pass"])

    def test_cli_exit_status(self):
        self.write_rows([self.row()])
        good = subprocess.run(
            [sys.executable, str(SCRIPT), str(self.csv)], capture_output=True, text=True
        )
        self.assertEqual(good.returncode, 0)
        self.assertTrue(json.loads(good.stdout)["format_and_disclosure_checks_pass"])
        self.write_rows([self.row(sequence="X")])
        bad = subprocess.run(
            [sys.executable, str(SCRIPT), str(self.csv)], capture_output=True, text=True
        )
        self.assertEqual(bad.returncode, 2)
        self.assertFalse(json.loads(bad.stdout)["format_and_disclosure_checks_pass"])


if __name__ == "__main__":
    unittest.main()
