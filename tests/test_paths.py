import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from fq_manifestor import fq_manifestor


class TestManifestPaths(TestCase):
    def setUp(self):
        temp_dir = TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        self.root = Path(temp_dir.name).resolve()
        self.addCleanup(os.chdir, os.getcwd())
        self.task = self.root / "first task"
        self.task.mkdir()
        os.chdir(self.task)

    def make_reads(self, directory):
        directory.mkdir(parents=True, exist_ok=True)
        for direction in ["R1", "R2"]:
            (directory / f"sampleA_{direction}_001.fastq.gz").touch()

    def manifest_paths(self, manifest):
        return manifest.read_text().splitlines()[1].split("\t")[1:]

    def test_default_paths_work_from_another_directory(self):
        reads = self.task / "reads"
        self.make_reads(reads)
        fq_manifestor("reads", "manifest.tsv", verbose=False)

        os.chdir(self.root)
        paths = self.manifest_paths(self.task / "manifest.tsv")
        self.assertEqual(paths, [
            str(reads / "sampleA_R1_001.fastq.gz"),
            str(reads / "sampleA_R2_001.fastq.gz"),
        ])
        for path in paths:
            self.assertTrue(Path(path).is_file())

    def test_paths_survive_relocation_to_another_task(self):
        self.make_reads(self.task / "reads" / "nested directory")
        fq_manifestor("reads", "manifest.tsv", verbose=False, use_pwd=True)

        second_task = self.root / "second task"
        shutil.copytree(self.task, second_task)
        os.chdir(second_task)
        shutil.rmtree(self.task)

        paths = self.manifest_paths(second_task / "manifest.tsv")
        self.assertEqual(paths, [
            "$PWD/reads/nested directory/sampleA_R1_001.fastq.gz",
            "$PWD/reads/nested directory/sampleA_R2_001.fastq.gz",
        ])
        for path in paths:
            # Resolve against the importing task, after the generating task is gone.
            self.assertTrue(Path(path.replace("$PWD", str(second_task))).is_file())

    def test_absolute_input_is_relative_to_cwd_not_manifest_directory(self):
        reads = self.task / "reads"
        self.make_reads(reads)
        (self.task / "results").mkdir()
        fq_manifestor(reads, "results/manifest.tsv", verbose=False, use_pwd=True)

        self.assertEqual(self.manifest_paths(Path("results/manifest.tsv")), [
            "$PWD/reads/sampleA_R1_001.fastq.gz",
            "$PWD/reads/sampleA_R2_001.fastq.gz",
        ])

    def test_staged_symlink_name_is_preserved(self):
        source = self.root / "original reads"
        self.make_reads(source)
        (self.task / "reads").symlink_to(source, target_is_directory=True)
        fq_manifestor("reads", "manifest.tsv", verbose=False, use_pwd=True)

        self.assertEqual(self.manifest_paths(Path("manifest.tsv")), [
            "$PWD/reads/sampleA_R1_001.fastq.gz",
            "$PWD/reads/sampleA_R2_001.fastq.gz",
        ])

    def test_input_outside_cwd_keeps_parent_components(self):
        self.make_reads(self.root / "reads")
        fq_manifestor("../reads", "manifest.tsv", verbose=False, use_pwd=True)

        self.assertEqual(self.manifest_paths(Path("manifest.tsv")), [
            "$PWD/../reads/sampleA_R1_001.fastq.gz",
            "$PWD/../reads/sampleA_R2_001.fastq.gz",
        ])

    def test_use_pwd_also_formats_unmatched_report(self):
        self.make_reads(self.task / "reads")
        for direction in ["R1", "R2"]:
            Path(f"reads/sampleB_{direction}_001.fastq.gz").touch()
        Path("metadata.tsv").write_text("sample-id\nsampleA\n")

        fq_manifestor("reads", "manifest.tsv", metadata="metadata.tsv",
                      verbose=False, use_pwd=True)

        self.assertEqual(self.manifest_paths(Path("manifest.tsv")), [
            "$PWD/reads/sampleA_R1_001.fastq.gz",
            "$PWD/reads/sampleA_R2_001.fastq.gz",
        ])
        self.assertEqual(Path("files-missing-from-md.tsv").read_text(),
                         "absolute-filepath\n"
                         "$PWD/reads/sampleB_R1_001.fastq.gz\n"
                         "$PWD/reads/sampleB_R2_001.fastq.gz\n")
