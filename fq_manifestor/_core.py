#!/usr/bin/env python

import csv
import glob
import os.path
import re


_METADATA_SAMPLE_ID_HEADERS = {'id', 'sampleid', 'sample id', 'sample-id'}
_MISSING_FILES_FILENAME = 'files-missing-from-md.tsv'


def _pwd_relative_path(filepath):
    """Defer the working-directory prefix until QIIME reads the manifest."""
    return os.path.join('$PWD', os.path.relpath(filepath))


def _read_metadata_sample_ids(metadata: str):
    '''
    Parses and returns all sample IDs from a metadata file.

    Parameters
    ----------
    metadata : str
        Path to the metadata file.

    Returns
    -------
    set[str]
        The parsed sample IDs.

    Raises
    ------
    ValueError
        If the metadata file is empty, no valid sample ID column names are
        found, or if more than one valid sample ID column name is found.
    '''
    with open(metadata, encoding='utf-8-sig', newline='') as md:
        reader = csv.reader(md, delimiter='\t')
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f'Metadata file is empty: {metadata}')

        normalized_header = [column.strip().lower() for column in header]
        id_columns = [
            index for index, column in enumerate(normalized_header)
            if column in _METADATA_SAMPLE_ID_HEADERS
        ]

        if len(id_columns) == 0:
            raise ValueError(
                'Metadata file does not contain a supported sample ID column. '
                'Expected one of: id, sampleid, sample id, sample-id.'
            )
        if len(id_columns) > 1:
            raise ValueError(
                'Metadata file contains multiple supported sample ID columns: '
                f'{(normalized_header[i] for i in id_columns)}.'
            )

        id_column = id_columns[0]
        sample_ids = set()
        for row in reader:
            if id_column >= len(row):
                continue

            sample_id = row[id_column].strip()
            if sample_id:
                sample_ids.add(sample_id)

    return sample_ids


def _check_metadata_overlap(fq_filepaths, metadata: str):
    '''
    Match FASTQ filenames to metadata IDs and validate file overlap.

    Parameters
    ----------
    fq_filepaths : list[str]
        All distinct FASTQ paths in the selected search scope.
    metadata : str
        Path to the metadata file.

    Returns
    -------
    list[tuple[str, str]]
        Matched sample IDs and FASTQ paths, in input file order.
    list[str]
        Unmatched FASTQ paths, in input file order.

    Raises
    ------
    ValueError
        If the metadata cannot be parsed, a filename contains multiple IDs,
        or more files are unmatched than matched. IDs are matched as literal,
        case-sensitive substrings of the filename, never the parent path.
    '''
    metadata_sample_ids = sorted(_read_metadata_sample_ids(metadata))
    matched_files = []
    unmatched_files = []
    for fq_filepath in fq_filepaths:
        fq_filename = os.path.basename(fq_filepath)
        matches = [sid for sid in metadata_sample_ids if sid in fq_filename]
        if len(matches) > 1:
            raise ValueError(
                f'Multiple metadata sample IDs match file {fq_filepath}: '
                f'{matches}. No output files were written.'
            )
        if matches:
            matched_files.append((matches[0], fq_filepath))
        else:
            unmatched_files.append(fq_filepath)

    if len(unmatched_files) > len(matched_files):
        raise ValueError(
            'There were more FASTQ files not represented in the metadata '
            f'({len(unmatched_files)}) than files represented in the metadata '
            f'({len(matched_files)}). No output files were written.'
        )

    return matched_files, unmatched_files


def fq_manifestor(input_dir,
                  output_fp,
                  fq_extensions=['fastq.gz', 'fq.gz'],
                  sample_id_pattern='(.*?)_',
                  f_read_pattern='_R1_',
                  r_read_pattern='_R2_',
                  filter_pattern=None,
                  metadata=None,
                  verbose=True,
                  use_pwd=False):
    """Write absolute read paths in a paired-end manifest and unmatched report.

    With use_pwd=True, write paths relative to the current directory, prefixed
    by literal $PWD. QIIME expands it at import time, so import from a directory
    with the same relative layout of read files.
    """

    input_dir = os.path.abspath(input_dir)
    if verbose:
        print("Searching directory: %s" % input_dir)

    fq_filepaths = []
    for fq_extension in fq_extensions:
        fq_filepaths += glob.glob(
            f'{input_dir}/**/*.{fq_extension}', recursive = True
        )

    if filter_pattern is not None:
        fq_filepaths = [fp for fp in fq_filepaths if filter_pattern in fp]

    # Count each real file once, including when extension patterns overlap.
    fq_filepaths = sorted({fp for fp in fq_filepaths if os.path.isfile(fp)})
    n_fq_filepaths = len(fq_filepaths)
    if verbose:
        print('Found %d fastq files.' % n_fq_filepaths)

    unmatched_files = None
    missing_files_fp = None
    if metadata is not None:
        missing_files_fp = os.path.join(
            os.path.dirname(os.path.abspath(output_fp)), _MISSING_FILES_FILENAME,
        )
        if len({os.path.realpath(fp) for fp in
                (metadata, output_fp, missing_files_fp)}) != 3:
            raise ValueError(
                'Metadata, manifest, and missing-files report paths must be '
                'distinct. No output files were written.'
            )
        files_with_ids, unmatched_files = _check_metadata_overlap(
            fq_filepaths, metadata,
        )
        if verbose:
            print(
                'Metadata check matched %d FASTQ files and excluded %d.' %
                (len(files_with_ids), len(unmatched_files))
            )
    else:
        files_with_ids = []
        for fq_filepath in fq_filepaths:
            fq_filename = os.path.basename(fq_filepath)
            m = re.search(sample_id_pattern, fq_filename)
            if not m or not m.group(1):
                raise ValueError('Sample ID not found in file: %s' % fq_filepath)
            files_with_ids.append((m.group(1), fq_filepath))

    sids_to_fps = {}
    for sid, fq_filepath in files_with_ids:
        fq_filename = os.path.basename(fq_filepath)
        if bool(re.search(f_read_pattern, fq_filename)):
            forward = True
            reverse = False
        elif bool(re.search(r_read_pattern, fq_filename)):
            forward = False
            reverse = True
        else:
            raise ValueError('Forward/reverse patterns not found in file: %s' % fq_filepath)

        try:
            sid_fps = sids_to_fps[sid]
        except KeyError:
            sid_fps = [None, None]

        read_index = 0 if forward else 1
        if metadata is not None and sid_fps[read_index] is not None:
            direction = 'forward' if forward else 'reverse'
            raise ValueError(
                f'Multiple {direction} reads for sample {sid}: '
                f'{sid_fps[read_index]}, {fq_filepath}. '
                'No output files were written.'
            )

        if forward:
            sid_fps[0] = fq_filepath
        else:
            sid_fps[1] = fq_filepath
        sids_to_fps[sid] = sid_fps

    records = []
    for sid, (fwd_fq_filepath, rev_fq_filepath) in sids_to_fps.items():
        if fwd_fq_filepath is None:
            raise ValueError('Missing forward read for sample: %s' % sid)

        if rev_fq_filepath is None:
            raise ValueError('Missing reverse read for sample: %s' % sid)

        records.append((sid, fwd_fq_filepath, rev_fq_filepath))

    if len(records) != (len(files_with_ids) / 2) and verbose:
        print("\n** WARNING**: "
              "The number of manifest records doesn't align with the number of "
              "fastq files that were found. It's possible that the match "
              "patterns aren't working correctly for your files. These can "
              "be customized when using the API.\n")

    format_path = _pwd_relative_path if use_pwd else os.path.abspath
    lines = ['sample-id\tforward-absolute-filepath\treverse-absolute-filepath']
    lines.extend(
        '%s\t%s\t%s' % (sid, format_path(fwd), format_path(rev))
        for sid, fwd, rev in records
    )

    with open(output_fp, 'w', encoding='utf-8') as of:
        of.write('\n'.join(lines))
        of.write('\n')

    if unmatched_files is not None:
        with open(missing_files_fp, 'w', encoding='utf-8', newline='') as of:
            writer = csv.writer(of, delimiter='\t', lineterminator='\n')
            writer.writerow(['absolute-filepath'])
            writer.writerows([format_path(fp)] for fp in unmatched_files)
