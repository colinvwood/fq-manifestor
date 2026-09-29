import os
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

import pytest

from fq_manifestor import fq_manifestor


def _touch(path):
    """Create an empty file at path, creating parent dirs as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()


def _read_manifest(path):
    lines = path.read_text().splitlines()
    assert lines[0] == "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"
    rows = {}
    for line in lines[1:]:
        if line:
            sid, fwd, rev = line.split("\t")
            rows[sid] = (fwd, rev)
    return rows


def _create_sample(tmp_path, sample_id):
    _touch(tmp_path / f"{sample_id}_S1_L001_R1_001.fastq.gz")
    _touch(tmp_path / f"{sample_id}_S1_L001_R2_001.fastq.gz")


def _write_metadata(path, header, sample_ids):
    path.write_text(
        "\n".join([f"{header}\tgroup", *[f"{sid}\ta" for sid in sample_ids]])
        + "\n"
    )


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

def test_basic_paired_end(tmp_path):
    _touch(tmp_path / "sampleA_S1_L001_R1_001.fastq.gz")
    _touch(tmp_path / "sampleA_S1_L001_R2_001.fastq.gz")
    _touch(tmp_path / "sampleB_S2_L001_R1_001.fastq.gz")
    _touch(tmp_path / "sampleB_S2_L001_R2_001.fastq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(str(tmp_path), str(out), verbose=False)

    rows = _read_manifest(out)
    assert set(rows) == {"sampleA", "sampleB"}
    assert os.path.isabs(rows["sampleA"][0])
    assert rows["sampleA"][0].endswith("sampleA_S1_L001_R1_001.fastq.gz")
    assert os.path.isabs(rows["sampleA"][1])
    assert rows["sampleA"][1].endswith("sampleA_S1_L001_R2_001.fastq.gz")


def test_recursive_search(tmp_path):
    sub = tmp_path / "sampleA_L001-ds.abc123"
    _touch(sub / "sampleA_S1_L001_R1_001.fastq.gz")
    _touch(sub / "sampleA_S1_L001_R2_001.fastq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(str(tmp_path), str(out), verbose=False)

    rows = _read_manifest(out)
    assert "sampleA" in rows


def test_filter_pattern(tmp_path):
    _touch(tmp_path / "keep_S1_L001_R1_001.fastq.gz")
    _touch(tmp_path / "keep_S1_L001_R2_001.fastq.gz")
    _touch(tmp_path / "drop_S2_L001_R1_001.fastq.gz")
    _touch(tmp_path / "drop_S2_L001_R2_001.fastq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(str(tmp_path), str(out), filter_pattern="keep", verbose=False)

    rows = _read_manifest(out)
    assert set(rows) == {"keep"}


def test_custom_extensions(tmp_path):
    _touch(tmp_path / "sampleA_S1_L001_R1_001.fq.gz")
    _touch(tmp_path / "sampleA_S1_L001_R2_001.fq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(str(tmp_path), str(out), fq_extensions=["fq.gz"], verbose=False)

    rows = _read_manifest(out)
    assert "sampleA" in rows


def test_custom_sample_id_pattern(tmp_path):
    _touch(tmp_path / "sampleA-S1-L001-R1-001.fastq.gz")
    _touch(tmp_path / "sampleA-S1-L001-R2-001.fastq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(
        str(tmp_path), str(out), sample_id_pattern="(.*?)-",
        f_read_pattern="-R1-", r_read_pattern="-R2-", verbose=False
    )

    rows = _read_manifest(out)
    assert "sampleA" in rows


def test_custom_read_patterns(tmp_path):
    _touch(tmp_path / "sampleA_S1_L001_1.fastq.gz")
    _touch(tmp_path / "sampleA_S1_L001_2.fastq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(
        str(tmp_path), str(out),
        f_read_pattern=r"_1\.fastq", r_read_pattern=r"_2\.fastq",
        verbose=False,
    )

    rows = _read_manifest(out)
    assert "sampleA" in rows


def test_verbose_false(tmp_path, capsys):
    _touch(tmp_path / "sampleA_S1_R1_001.fastq.gz")
    _touch(tmp_path / "sampleA_S1_R2_001.fastq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(str(tmp_path), str(out), verbose=False)

    captured = capsys.readouterr()
    assert captured.out == ""


def test_output_file_written(tmp_path):
    _touch(tmp_path / "sampleA_S1_R1_001.fastq.gz")
    _touch(tmp_path / "sampleA_S1_R2_001.fastq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(str(tmp_path), str(out), verbose=False)

    assert out.exists()
    assert out.stat().st_size > 0


# ---------------------------------------------------------------------------
# Metadata checks
# ---------------------------------------------------------------------------

class TestMetadataParameter(TestCase):
    def setUp(self):
        temp_dir = TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        self.tmp_path = Path(temp_dir.name)

    def assertManifestSampleIds(self, path, expected_ids):
        lines = path.read_text().splitlines()
        self.assertEqual(
            lines[0],
            "sample-id\tforward-absolute-filepath\treverse-absolute-filepath",
        )
        sample_ids = set()
        for line in lines[1:]:
            if line:
                sid, fwd, rev = line.split("\t")
                sample_ids.add(sid)
        self.assertEqual(sample_ids, expected_ids)

    def test_metadata_ids_replace_sample_id_regex(self):
        _create_sample(self.tmp_path, "sample_with_underscores")
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sample_with_underscores"])
        out = self.tmp_path / "manifest.tsv"

        fq_manifestor(
            str(self.tmp_path), str(out), metadata=str(metadata),
            sample_id_pattern="[", verbose=False,
        )

        self.assertManifestSampleIds(out, {"sample_with_underscores"})

    def test_metadata_threshold_counts_unidentified_files(self):
        _create_sample(self.tmp_path, "sampleA")
        for name in ["unidentified.fastq.gz", "unpaired.fastq.gz",
                     "unrecognized.fastq.gz"]:
            _touch(self.tmp_path / name)
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sampleA"])
        out = self.tmp_path / "manifest.tsv"

        with self.assertRaisesRegex(
                ValueError, r"files not represented in the metadata \(3\).*"
                r"files represented in the metadata \(2\)"):
            fq_manifestor(
                str(self.tmp_path), str(out), metadata=str(metadata),
                verbose=False,
            )

        self.assertFalse(out.exists())

    def test_metadata_filters_manifest_and_writes_missing_files(self):
        for index, sample_id_header in enumerate([
                "id", "sampleid", "sample id", "sample-id",
                "SampleID", "Sample ID"]):
            with self.subTest(sample_id_header=sample_id_header):
                tmp_path = self.tmp_path / str(index)
                for sample_id in ["sampleA", "sampleB", "sampleC"]:
                    _create_sample(tmp_path, sample_id)

                metadata = tmp_path / "metadata.tsv"
                _write_metadata(
                    metadata, sample_id_header, ["sampleB", "sampleC"]
                )
                out = tmp_path / "manifest.tsv"

                fq_manifestor(
                    str(tmp_path), str(out), metadata=str(metadata),
                    verbose=False,
                )

                self.assertManifestSampleIds(out, {"sampleB", "sampleC"})
                missing = tmp_path / "files-missing-from-md.tsv"
                self.assertEqual(
                    missing.read_text().splitlines(),
                    ["absolute-filepath",
                     str(tmp_path / "sampleA_S1_L001_R1_001.fastq.gz"),
                     str(tmp_path / "sampleA_S1_L001_R2_001.fastq.gz")],
                )

    def test_metadata_all_samples_present_writes_empty_report(self):
        _create_sample(self.tmp_path, "sampleA")
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sampleA", "sampleA", "extra"])
        out = self.tmp_path / "manifest.tsv"

        fq_manifestor(
            str(self.tmp_path), str(out), metadata=str(metadata), verbose=False
        )

        self.assertManifestSampleIds(out, {"sampleA"})
        missing = self.tmp_path / "files-missing-from-md.tsv"
        self.assertEqual(missing.read_text(), "absolute-filepath\n")

    def test_metadata_allows_equal_removed_and_retained(self):
        _create_sample(self.tmp_path, "sampleA")
        _create_sample(self.tmp_path, "sampleB")
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "id", ["sampleA"])
        out = self.tmp_path / "manifest.tsv"

        fq_manifestor(
            str(self.tmp_path), str(out), metadata=str(metadata), verbose=False
        )

        self.assertManifestSampleIds(out, {"sampleA"})
        missing = self.tmp_path / "files-missing-from-md.tsv"
        self.assertEqual(
            missing.read_text().splitlines(),
            ["absolute-filepath",
             str(self.tmp_path / "sampleB_S1_L001_R1_001.fastq.gz"),
             str(self.tmp_path / "sampleB_S1_L001_R2_001.fastq.gz")],
        )

    def test_metadata_rejects_removing_more_files_than_retained(self):
        for sample_id in ["sampleA", "sampleB", "sampleC"]:
            _create_sample(self.tmp_path, sample_id)

        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sampleA"])
        out = self.tmp_path / "manifest.tsv"

        with self.assertRaisesRegex(
                ValueError,
                r"not represented in the metadata \(4\) than files represented "
                r"in the metadata \(2\)"):
            fq_manifestor(
                str(self.tmp_path), str(out), metadata=str(metadata),
                verbose=False,
            )

        self.assertFalse(out.exists())
        self.assertFalse(
            (self.tmp_path / "files-missing-from-md.tsv").exists()
        )

    def test_metadata_rejects_unsupported_sample_id_header(self):
        for header in ["#SampleID", "feature-id", "subject"]:
            with self.subTest(header=header):
                tmp_path = self.tmp_path / header
                _create_sample(tmp_path, "sampleA")
                metadata = tmp_path / "metadata.tsv"
                _write_metadata(metadata, header, ["sampleA"])
                out = tmp_path / "manifest.tsv"

                with self.assertRaisesRegex(
                        ValueError, "does not contain a supported"):
                    fq_manifestor(
                        str(tmp_path), str(out), metadata=str(metadata),
                        verbose=False,
                    )

                self.assertFalse(out.exists())

    def test_metadata_rejects_multiple_sample_id_columns(self):
        _create_sample(self.tmp_path, "sampleA")
        metadata = self.tmp_path / "metadata.tsv"
        metadata.write_text("id\tsample-id\nsampleA\tsampleA\n")
        out = self.tmp_path / "manifest.tsv"

        with self.assertRaisesRegex(ValueError, "multiple supported"):
            fq_manifestor(
                str(self.tmp_path), str(out), metadata=str(metadata),
                verbose=False,
            )

    def test_metadata_sample_id_matching_is_case_sensitive(self):
        _create_sample(self.tmp_path, "sampleA")
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["SampleA"])
        out = self.tmp_path / "manifest.tsv"

        with self.assertRaisesRegex(
                ValueError,
                r"not represented in the metadata \(2\) than files represented "
                r"in the metadata \(0\)"):
            fq_manifestor(
                str(self.tmp_path), str(out), metadata=str(metadata),
                verbose=False,
            )

    def test_metadata_matches_literal_id_anywhere_in_filename(self):
        sample_id = "sample_[A]+.42"
        paths = [self.tmp_path / f"run-prefix-{sample_id}_R{read}_001.fastq.gz"
                 for read in [1, 2]]
        for path in paths:
            _touch(path)
        # This would match if the ID were interpreted as a regex.
        unrelated = self.tmp_path / "sample_Ax42.fastq.gz"
        _touch(unrelated)
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", [sample_id])
        out = self.tmp_path / "manifest.tsv"

        fq_manifestor(
            str(self.tmp_path), str(out), metadata=str(metadata), verbose=False,
        )

        self.assertEqual(_read_manifest(out), {sample_id: tuple(map(str, paths))})
        self.assertEqual(
            (self.tmp_path / "files-missing-from-md.tsv").read_text(),
            f"absolute-filepath\n{unrelated}\n",
        )

    def test_metadata_rejects_ambiguous_id_matches_before_writing(self):
        _create_sample(self.tmp_path, "sample10")
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sample1", "sample10"])
        out = self.tmp_path / "manifest.tsv"

        with self.assertRaisesRegex(ValueError, "Multiple metadata sample IDs"):
            fq_manifestor(
                str(self.tmp_path), str(out), metadata=str(metadata),
                verbose=False,
            )

        self.assertFalse(out.exists())
        self.assertFalse((self.tmp_path / "files-missing-from-md.tsv").exists())

    def test_metadata_does_not_match_parent_directory(self):
        _create_sample(self.tmp_path / "sampleA", "unrelated")
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sampleA"])
        out = self.tmp_path / "manifest.tsv"

        with self.assertRaisesRegex(
                ValueError, r"files represented in the metadata \(0\)"):
            fq_manifestor(
                str(self.tmp_path), str(out), metadata=str(metadata),
                verbose=False,
            )

        self.assertFalse(out.exists())

    def test_metadata_allows_half_unidentified_files(self):
        _create_sample(self.tmp_path, "sampleA")
        missing_paths = [self.tmp_path / "unidentified.fastq.gz",
                         self.tmp_path / "unpaired.fastq.gz"]
        for path in missing_paths:
            _touch(path)
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sampleA"])
        out = self.tmp_path / "manifest.tsv"

        fq_manifestor(
            str(self.tmp_path), str(out), metadata=str(metadata), verbose=False,
        )

        self.assertManifestSampleIds(out, {"sampleA"})
        self.assertEqual(
            (self.tmp_path / "files-missing-from-md.tsv").read_text().splitlines(),
            ["absolute-filepath", *map(str, missing_paths)],
        )

    def test_metadata_filter_pattern_limits_search_scope(self):
        _create_sample(self.tmp_path / "selected", "sampleA")
        for sample_id in ["sampleB", "sampleC", "sampleD"]:
            _create_sample(self.tmp_path / "excluded", sample_id)
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sampleA"])
        out = self.tmp_path / "manifest.tsv"

        fq_manifestor(
            str(self.tmp_path), str(out), metadata=str(metadata),
            filter_pattern="selected", verbose=False,
        )

        self.assertManifestSampleIds(out, {"sampleA"})
        self.assertEqual(
            (self.tmp_path / "files-missing-from-md.tsv").read_text(),
            "absolute-filepath\n",
        )

    def test_metadata_counts_distinct_files_with_selected_extensions(self):
        _create_sample(self.tmp_path, "sampleA")
        unmatched = self.tmp_path / "unidentified.fq.gz"
        _touch(unmatched)
        (self.tmp_path / "not-a-file.fastq.gz").mkdir()
        _touch(self.tmp_path / "unselected.txt")
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sampleA"])
        out = self.tmp_path / "manifest.tsv"

        fq_manifestor(
            str(self.tmp_path), str(out), metadata=str(metadata),
            fq_extensions=["fastq.gz", "fq.gz", "gz"], verbose=False,
        )

        self.assertManifestSampleIds(out, {"sampleA"})
        self.assertEqual(
            (self.tmp_path / "files-missing-from-md.tsv").read_text(),
            f"absolute-filepath\n{unmatched}\n",
        )

    def test_metadata_still_validates_matched_read_pairs(self):
        for read, message in [("R1", "Missing reverse"),
                              ("R2", "Missing forward"),
                              ("R3", "Forward/reverse patterns not found")]:
            with self.subTest(read=read):
                root = self.tmp_path / read
                _touch(root / f"sampleA_{read}_001.fastq.gz")
                metadata = root / "metadata.tsv"
                _write_metadata(metadata, "sample-id", ["sampleA"])
                out = root / "manifest.tsv"

                with self.assertRaisesRegex(ValueError, message):
                    fq_manifestor(
                        str(root), str(out), metadata=str(metadata),
                        verbose=False,
                    )

                self.assertFalse(out.exists())
                self.assertFalse((root / "files-missing-from-md.tsv").exists())

    def test_metadata_rejects_duplicate_read_direction(self):
        for read in ["R1", "R2"]:
            with self.subTest(read=read):
                root = self.tmp_path / read
                _create_sample(root, "sampleA")
                _touch(root / f"sampleA_S2_L002_{read}_001.fastq.gz")
                metadata = root / "metadata.tsv"
                _write_metadata(metadata, "sample-id", ["sampleA"])
                out = root / "manifest.tsv"

                with self.assertRaisesRegex(ValueError, "Multiple .* reads"):
                    fq_manifestor(
                        str(root), str(out), metadata=str(metadata),
                        verbose=False,
                    )

                self.assertFalse(out.exists())
                self.assertFalse((root / "files-missing-from-md.tsv").exists())

    def test_metadata_threshold_failure_preserves_existing_outputs(self):
        _create_sample(self.tmp_path, "sampleA")
        for name in ["unidentified", "unpaired", "unrecognized"]:
            _touch(self.tmp_path / f"{name}.fastq.gz")
        metadata = self.tmp_path / "metadata.tsv"
        _write_metadata(metadata, "sample-id", ["sampleA"])
        out = self.tmp_path / "manifest.tsv"
        report = self.tmp_path / "files-missing-from-md.tsv"
        out.write_text("previous manifest\n")
        report.write_text("previous report\n")

        with self.assertRaisesRegex(ValueError, "more FASTQ files"):
            fq_manifestor(
                str(self.tmp_path), str(out), metadata=str(metadata),
                verbose=False,
            )

        self.assertEqual(out.read_text(), "previous manifest\n")
        self.assertEqual(report.read_text(), "previous report\n")

    def test_metadata_rejects_colliding_output_paths(self):
        _create_sample(self.tmp_path, "sampleA")
        for metadata_name, output_name in [
                ("metadata.tsv", "metadata.tsv"),
                ("metadata.tsv", "files-missing-from-md.tsv"),
                ("files-missing-from-md.tsv", "manifest.tsv")]:
            with self.subTest(metadata=metadata_name, output=output_name):
                metadata = self.tmp_path / metadata_name
                _write_metadata(metadata, "sample-id", ["sampleA"])
                original = metadata.read_text()

                with self.assertRaisesRegex(ValueError, "paths must be distinct"):
                    fq_manifestor(
                        str(self.tmp_path), str(self.tmp_path / output_name),
                        metadata=str(metadata), verbose=False,
                    )

                self.assertEqual(metadata.read_text(), original)


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

def test_missing_forward_read(tmp_path):
    _touch(tmp_path / "sampleA_S1_R2_001.fastq.gz")  # only reverse

    out = tmp_path / "manifest.tsv"
    with pytest.raises(ValueError, match="Missing forward read for sample: sampleA"):
        fq_manifestor(str(tmp_path), str(out), verbose=False)


def test_missing_reverse_read(tmp_path):
    _touch(tmp_path / "sampleA_S1_R1_001.fastq.gz")  # only forward

    out = tmp_path / "manifest.tsv"
    with pytest.raises(ValueError, match="Missing reverse read for sample: sampleA"):
        fq_manifestor(str(tmp_path), str(out), verbose=False)


def test_no_read_pattern_match(tmp_path):
    # File name matches split pattern but neither R1 nor R2 pattern
    _touch(tmp_path / "sampleA_S1_L001_R3_001.fastq.gz")

    out = tmp_path / "manifest.tsv"
    with pytest.raises(ValueError, match="Forward/reverse patterns not found"):
        fq_manifestor(str(tmp_path), str(out), verbose=False)


# ---------------------------------------------------------------------------
# Warning case
# ---------------------------------------------------------------------------

def test_warning_mismatch(tmp_path, capsys):
    # Three files: two R1 for the same sample ID (second overwrites first),
    # plus one R2. sids_to_fps ends up with 1 complete record but 3 files
    # were found, so (1 != 3/2) triggers the WARNING.
    _touch(tmp_path / "sampleA_S1_R1_001.fastq.gz")
    _touch(tmp_path / "sampleA_S2_R1_002.fastq.gz")  # same sample id, overwrites
    _touch(tmp_path / "sampleA_S1_R2_001.fastq.gz")

    out = tmp_path / "manifest.tsv"
    fq_manifestor(str(tmp_path), str(out), verbose=True)

    captured = capsys.readouterr()
    assert "WARNING" in captured.out
