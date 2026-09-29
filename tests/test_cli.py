import os
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from typer.testing import CliRunner

from fq_manifestor.cli import app


class TestCLI(TestCase):
    def setUp(self):
        temp_dir = TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        self.tmp_path = Path(temp_dir.name)
        self.runner = CliRunner()

    def test_default_paths_are_absolute(self):
        for direction in ["R1", "R2"]:
            (self.tmp_path / f"sampleA_{direction}_001.fastq.gz").touch()
        output = self.tmp_path / "manifest.tsv"

        result = self.runner.invoke(app, [str(self.tmp_path), str(output)])

        self.assertEqual(result.exit_code, 0, msg=str(result.exception))
        self.assertEqual(output.read_text().splitlines()[1],
                         f"sampleA\t{self.tmp_path}/sampleA_R1_001.fastq.gz\t"
                         f"{self.tmp_path}/sampleA_R2_001.fastq.gz")
        self.assertIn("To import, try running:", result.stdout)
        self.assertNotIn("same relative read paths", result.stdout)

    def test_use_pwd_option(self):
        self.addCleanup(os.chdir, os.getcwd())
        os.chdir(self.tmp_path)
        reads = Path("reads")
        reads.mkdir()
        for direction in ["R1", "R2"]:
            (reads / f"sampleA_{direction}_001.fastq.gz").touch()

        result = self.runner.invoke(app, ["reads", "manifest.tsv", "--use-pwd"])

        self.assertEqual(result.exit_code, 0, msg=str(result.exception))
        self.assertEqual(Path("manifest.tsv").read_text().splitlines()[1],
                         "sampleA\t$PWD/reads/sampleA_R1_001.fastq.gz\t"
                         "$PWD/reads/sampleA_R2_001.fastq.gz")
        self.assertIn("same relative read paths", result.stdout)

    def test_metadata_option(self):
        (self.tmp_path / "run_sample_A_S1_R1_001.fastq.gz").touch()
        (self.tmp_path / "run_sample_A_S1_R2_001.fastq.gz").touch()
        (self.tmp_path / "sampleB_S2_R1_001.fastq.gz").touch()
        (self.tmp_path / "sampleB_S2_R2_001.fastq.gz").touch()
        metadata = self.tmp_path / "metadata.tsv"
        metadata.write_text("sample-id\nsample_A\n")
        output = self.tmp_path / "manifest.tsv"

        result = self.runner.invoke(
            app,
            [str(self.tmp_path), str(output), "--metadata", str(metadata),
             "--sample-id-pattern", "[", "--no-verbose"],
        )

        self.assertEqual(result.exit_code, 0, msg=str(result.exception))
        self.assertTrue(output.read_text().splitlines()[1].startswith("sample_A\t"))
        self.assertNotIn("sampleB", output.read_text())
        missing = self.tmp_path / "files-missing-from-md.tsv"
        self.assertEqual(
            missing.read_text().splitlines(),
            ["absolute-filepath",
             str(self.tmp_path / "sampleB_S2_R1_001.fastq.gz"),
             str(self.tmp_path / "sampleB_S2_R2_001.fastq.gz")],
        )

    def test_metadata_threshold_error(self):
        for name in ["sampleA_R1_001", "sampleA_R2_001", "unknown1",
                     "unknown2", "unknown3"]:
            (self.tmp_path / f"{name}.fastq.gz").touch()
        metadata = self.tmp_path / "metadata.tsv"
        metadata.write_text("sample-id\nsampleA\n")
        output = self.tmp_path / "manifest.tsv"

        result = self.runner.invoke(
            app, [str(self.tmp_path), str(output), "--metadata", str(metadata)],
        )

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("more FASTQ files", str(result.exception))
        self.assertFalse(output.exists())
        self.assertFalse((self.tmp_path / "files-missing-from-md.tsv").exists())
