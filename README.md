# fq-manifestor

`fq-manifestor` is a tool to assist with generating [QIIME 2 paired-end fastq manifest files](https://amplicon-docs.qiime2.org/en/stable/how-to-guides/how-to-import.html#import-fastq-manifest). This is not widely tested and likely not very general purpose, but is designed as a starting point to assist Illumina users with importing paired-end demultiplexed fastq files for use with QIIME 2. This tool does not actually import your data into QIIME 2, but rather generates a fastq manifest file and suggests an import command to use with QIIME 2. That increases portability of the tool as it can be run on systems without QIIME 2 installed on them.

This tool was originally developed to facilitate work with the [NAU Genetics Core Facility](https://in.nau.edu/gcf/).

## Why is this tool needed?

Importing can be hard. I talk a little about why [here](https://amplicon-docs.qiime2.org/en/stable/explanations/why-importing.html).

## Installation

```
pip install https://github.com/gregcaporaso/fq-manifestor/archive/refs/heads/main.zip
```

For development:

```
git clone https://github.com/gregcaporaso/fq-manifestor
cd fq-manifestor
pip install -e ".[dev]"
```

Requires Python 3.10–3.12.
Should be easily installable in any rachis deployments (e.g., QIIME 2, MOSHPIT) ≥ 2026.1. 

## Input data format

If you have a directory of files, `/sequence-data/humanure-run1`, containing demultiplexed paired end Illumina sequence data this tool will assume the directory is structured as follows:

```
$ tree /sequence-data/humanure-run1
/sequence-data/humanure-run1
├── 07717ac5_L001-ds.0795510ce4044fbab87b5c0b08aec105
│   ├── 07717ac5_S71_L001_R1_001.fastq.gz
│   └── 07717ac5_S71_L001_R2_001.fastq.gz
├── 08a10d08_L001-ds.060074faec87417b88e83d07470e7b10
│   ├── 08a10d08_S30_L001_R1_001.fastq.gz
│   └── 08a10d08_S30_L001_R2_001.fastq.gz
├── 091cbf06_L001-ds.9e1f729478e645bfa87b695aec93300d
│   ├── 091cbf06_S35_L001_R1_001.fastq.gz
│   └── 091cbf06_S35_L001_R2_001.fastq.gz
├── 10583739_1_L001-ds.0e7c88a021984a38b930ce5c60ddadfd
│   ├── 10583739-1_S72_L001_R1_001.fastq.gz
│   └── 10583739-1_S72_L001_R2_001.fastq.gz
```

In this case there are four samples, `07717ac5`, `08a10d08`, `091cbf06`, and `10583739-1`. The samples are grouped into per-sample directories, and there is a forward (R1) and reverse (R2) read file for each sample. Without metadata, this tool will assume that the text before the first `_` in the fastq filenames (not the directory names) is the sample id. Note that the sample ID preceding the first `_` must be the same for both forward (R1) and reverse (R2) read. With `--metadata`, the metadata IDs are instead matched directly against filenames, as described below.

For example, the sample id that will be associated with the following files:

```
├── 10583739_1_L001-ds.0e7c88a021984a38b930ce5c60ddadfd
│   ├── 10583739-1_S72_L001_R1_001.fastq.gz
│   └── 10583739-1_S72_L001_R2_001.fastq.gz
```

is `10583739-1`. This is the text before the first `_` in the file name. Notice that the `-` in the sample id is replaced with an `_` in the directory name - this tool is building sample ids from the file names, not the directory names.

## Usage

```
fq-manifestor [OPTIONS] INPUT_DIR OUTPUT_FP
```

**Arguments:**

- `INPUT_DIR` — Directory to search recursively for fastq files.
- `OUTPUT_FP` — Path to write the manifest TSV file.

**Options:**

| Option | Default | Description |
| --- | --- | --- |
| `--fq-extensions TEXT` | `fastq.gz,fq.gz` | Comma-separated list of fastq file extensions to search for. |
| `--sample-id-pattern TEXT` | `(.*?)_` | Regex with a capture group used to extract the sample ID from filenames. Ignored when `--metadata` is provided. |
| `--f-read-pattern TEXT` | `_R1_` | Regex pattern identifying forward reads. |
| `--r-read-pattern TEXT` | `_R2_` | Regex pattern identifying reverse reads. |
| `--filter-pattern TEXT` | | Optional substring; only files containing this string are included. |
| `--metadata PATH` | | Optional metadata TSV whose sample IDs are matched literally against filenames instead of using `--sample-id-pattern`. |
| `--use-pwd` | Off | Write read paths with a literal `$PWD` prefix instead of absolute paths, in both the manifest and unmatched-files report. |
| `--verbose / --no-verbose` | `--verbose` | Print progress messages. |
| `--help` | | Show help and exit. |

**Example:**

```
fq-manifestor /sequence-data/humanure-run1 humanure-run1-fq-manifest.txt
```

This will either error with a hopefully-informative error message (but no guarantees!) or succeed with output that looks something like:

```
Searching directory: /sequence-data/humanure-run1
Found 160 fastq files.
To import, try running:
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path humanure-run1-fq-manifest.txt --output-path demux.qza --input-format PairedEndFastqManifestPhred33V2
```

Try running that command in a QIIME 2 conda environment, and it should import successfully.

### Paths in the manifest

Read paths are absolute by default, so importing does not depend on the current
working directory as long as the reads remain at those locations.

Pass `--use-pwd` to write literal `$PWD` followed by the path relative to the
working directory where `fq-manifestor` ran. For example, running
`fq-manifestor reads manifest.tsv --use-pwd` produces entries such as
`$PWD/reads/sampleA_R1_001.fastq.gz`. QIIME expands `$PWD` when importing;
fq-manifestor does not expand it when writing. The QIIME-required
`forward-absolute-filepath` and `reverse-absolute-filepath` headers stay the same.

With `--use-pwd`, run the import from the same working directory, or another with
the same relative read layout. The manifest's own location does not determine
where read paths resolve. In Nextflow, both tasks can stage the input directory
as `reads`, so the manifest resolves against the importing task's staged reads.

`INPUT_DIR` can be absolute or relative in either mode. With `--use-pwd`, reads
outside the working directory produce paths containing `..`, so that relative
layout must also be preserved when moving the data.

### Finding samples using metadata

Pass a metadata TSV with `--metadata` to search fastq filenames directly for
its sample IDs. In this mode, `--sample-id-pattern` is ignored:

```
fq-manifestor /sequence-data/humanure-run1 humanure-run1-fq-manifest.txt \
    --metadata sample-metadata.tsv
```

The metadata must contain exactly one sample ID column named `id`, `sampleid`,
`sample id`, or `sample-id`. Header matching is case-insensitive, while sample
ID matching is case-sensitive. Empty and duplicate metadata IDs are ignored.

Each ID is matched as literal text anywhere in the filename, not in parent
directory names. Underscores and regex metacharacters in IDs have no special
meaning. For example, the metadata ID `sample_A` matches
`run-prefix_sample_A_S1_L001_R1_001.fastq.gz` and its corresponding R2 file,
and the manifest uses `sample_A` as the sample ID.

A file matching multiple metadata IDs causes an error rather than choosing an
ID. For example, a filename containing `sample10` is ambiguous if the metadata
contains both `sample1` and `sample10`. For matched files, `--f-read-pattern`
and `--r-read-pattern` still identify read direction. Missing mates, unidentified
read directions, and multiple files for the same sample and direction cause an
error. Metadata IDs without matching fastq files are ignored.

The comparison set consists of all distinct files matching `--fq-extensions`,
within the search scope selected by `--filter-pattern` when supplied. A file
that matches no metadata ID is excluded, even if its sample ID or read direction
cannot be inferred. If more files are unmatched than matched, the command raises
an error before writing either output. Exactly half unmatched is allowed. For
example, two matched read files and two unmatched files are allowed; two matched
read files and three unmatched files cause an error. This counts files, not
samples or manifest rows.

Unmatched file paths are written to `files-missing-from-md.tsv` beside the
manifest, with the column header `absolute-filepath`. They are absolute by
default and use literal `$PWD` when `--use-pwd` is supplied. This report replaces
`sample-ids-missing-from-md.tsv`, because unmatched sample IDs are no longer
inferred from filenames. The report contains only its header when every file
matches. The metadata, manifest, and report must have distinct paths.

## Using the Python API

```python
from fq_manifestor import fq_manifestor

fq_manifestor(
    input_dir='/sequence-data/humanure-run1',
    output_fp='humanure-run1-fq-manifest.txt',
    f_read_pattern='_R1_',
    r_read_pattern='_R2_',
    metadata='sample-metadata.tsv',
)
```

All parameters available through the CLI are also available in the Python API. See the source for full details.

Set `use_pwd=True` to enable `$PWD` paths in the Python API; the default is `False`.

## What if it doesn't work?

Let me know on the [issue tracker](https://github.com/gregcaporaso/fq-manifestor/issues), and I may be able to help. But this is very experimental so no promises!
