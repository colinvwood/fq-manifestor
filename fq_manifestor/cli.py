from __future__ import annotations

import os.path
from typing import Annotated, Optional

import typer

from fq_manifestor._core import fq_manifestor

app = typer.Typer(
    help="Generate a QIIME 2 paired-end fastq manifest from a directory of demultiplexed reads."
)


@app.command()
def main(
    input_dir: Annotated[
        str,
        typer.Argument(help="Directory to search recursively for fastq files."),
    ],
    output_fp: Annotated[
        str,
        typer.Argument(help="Path to write the manifest TSV file."),
    ],
    fq_extensions: Annotated[
        str,
        typer.Option(
            "--fq-extensions",
            help="Comma-separated list of fastq file extensions to search for.",
        ),
    ] = "fastq.gz,fq.gz",
    sample_id_pattern: Annotated[
        str,
        typer.Option(
            "--sample-id-pattern",
            help=("Regex with a capture group used to extract the sample ID "
                  "from filenames. Ignored when --metadata is provided."),
        ),
    ] = "(.*?)_",
    f_read_pattern: Annotated[
        str,
        typer.Option(
            "--f-read-pattern",
            help="Regex pattern identifying forward reads.",
        ),
    ] = "_R1_",
    r_read_pattern: Annotated[
        str,
        typer.Option(
            "--r-read-pattern",
            help="Regex pattern identifying reverse reads.",
        ),
    ] = "_R2_",
    filter_pattern: Annotated[
        Optional[str],
        typer.Option(
            "--filter-pattern",
            help="Optional substring; only files containing this string are included.",
        ),
    ] = None,
    metadata: Annotated[
        Optional[str],
        typer.Option(
            "--metadata",
            help=("Optional metadata TSV whose sample IDs are matched literally "
                  "against filenames instead of using --sample-id-pattern."),
        ),
    ] = None,
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose/--no-verbose",
            help="Print progress messages.",
        ),
    ] = True,
    use_pwd: Annotated[
        bool,
        typer.Option(
            "--use-pwd",
            help=("Write read paths using literal $PWD plus paths relative to "
                  "the current directory, instead of absolute paths. Applies "
                  "to the manifest and unmatched-files report."),
        ),
    ] = False,
) -> None:
    extensions = [e.strip() for e in fq_extensions.split(",") if e.strip()]

    fq_manifestor(
        input_dir=input_dir,
        output_fp=output_fp,
        fq_extensions=extensions,
        sample_id_pattern=sample_id_pattern,
        f_read_pattern=f_read_pattern,
        r_read_pattern=r_read_pattern,
        filter_pattern=filter_pattern,
        metadata=metadata,
        verbose=verbose,
        use_pwd=use_pwd,
    )

    abs_output_fp = os.path.abspath(output_fp)
    if use_pwd:
        typer.echo(
            "To import, run from this directory or one with the same relative "
            "read paths:"
        )
    else:
        typer.echo("To import, try running:")
    typer.echo(
        f"qiime tools import "
        f"--type 'SampleData[PairedEndSequencesWithQuality]' "
        f"--input-path {abs_output_fp} --output-path demux.qza "
        f"--input-format PairedEndFastqManifestPhred33V2"
    )
