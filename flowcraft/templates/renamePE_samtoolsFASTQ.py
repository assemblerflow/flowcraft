#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import sys
import bz2
import gzip
import zipfile
from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

"""
Purpose
-------

This module renames the fastq headers with PE terminations
that were not include in samtools fastq command


Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id``: Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fastq_pair`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'`

Generated output
----------------
- ``fastq_pair`` : Pair of FastQ file paths with rename headers.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'`

Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "09.09.2019"
__template__ = "retrieved_mapped-nf"

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    FASTQ_PAIR = '$fastq_pair'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))


COPEN = {
    "gz": gzip.open,
    "bz2": bz2.open,
    "zip": zipfile.ZipFile
}

MAGIC_DICT = {
    b"\\x1f\\x8b\\x08": "gz",
    b"\\x42\\x5a\\x68": "bz2",
    b"\\x50\\x4b\\x03\\x04": "zip"
}


def guess_file_compression(file_path, magic_dict=None):
    """Guesses the compression of an input file.

    This function guesses the compression of a given file by checking for
    a binary signature at the beginning of the file. These signatures are
    stored in the :py:data:`MAGIC_DICT` dictionary. The supported compression
    formats are gzip, bzip2 and zip. If none of the signatures in this
    dictionary are found at the beginning of the file, it returns ``None``.

    Parameters
    ----------
    file_path : str
        Path to input file.
    magic_dict : dict, optional
        Dictionary containing the signatures of the compression types. The
        key should be the binary signature and the value should be the
        compression format. If left ``None``, it falls back to
        :py:data:`MAGIC_DICT`.

    Returns
    -------
    file_type : str or None
        If a compression type is detected, returns a string with the format.
        If not, returns ``None``.
    """

    if not magic_dict:
        magic_dict = MAGIC_DICT

    max_len = max(len(x) for x in magic_dict)

    with open(file_path, "rb") as f:
        file_start = f.read(max_len)

    logger.debug("Binary signature start: {}".format(file_start))

    for magic, file_type in magic_dict.items():
        if file_start.startswith(magic):
            return file_type

    return None


def formartFastqHeaders(sample_name, in_fastq_1, in_fastq_2):
    out_fastq_1 = os.path.join(os.getcwd(), sample_name + '.headersRenamed_1.fq')
    out_fastq_2 = os.path.join(os.getcwd(), sample_name + '.headersRenamed_2.fq')

    writer_in_fastq_1 = open(out_fastq_1, 'wt')
    writer_in_fastq_2 = open(out_fastq_2, 'wt')

    outfiles = [out_fastq_1, out_fastq_2]

    with open(in_fastq_1, 'r') as reader_in_fastq_1, open(in_fastq_2, 'r') as reader_in_fastq_2:
        plus_line = True
        quality_line = True
        number_reads = 0
        for in_1, in_2 in zip(reader_in_fastq_1, reader_in_fastq_2):
            if len(in_1) > 0:
                in_1 = in_1.splitlines()[0]
                in_2 = in_2.splitlines()[0]

                if in_1.startswith('@') and plus_line and quality_line:
                    if in_1 != in_2:
                        sys.exit('The PE fastq files are not aligned properly!')
                    in_1 += '/1' + '\\n'
                    in_2 += '/2' + '\\n'
                    writer_in_fastq_1.write(in_1)
                    writer_in_fastq_2.write(in_2)
                    plus_line = False
                    quality_line = False
                elif in_1.startswith('+') and not plus_line:
                    in_1 += '\\n'
                    writer_in_fastq_1.write(in_1)
                    writer_in_fastq_2.write(in_1)
                    plus_line = True
                elif plus_line and not quality_line:
                    in_1 += '\\n'
                    in_2 += '\\n'
                    writer_in_fastq_1.write(in_1)
                    writer_in_fastq_2.write(in_2)
                    writer_in_fastq_1.flush()
                    writer_in_fastq_2.flush()
                    number_reads += 1
                    quality_line = True
                else:
                    in_1 += '\\n'
                    in_2 += '\\n'
                    writer_in_fastq_1.write(in_1)
                    writer_in_fastq_2.write(in_2)

    writer_in_fastq_1.close()
    writer_in_fastq_2.close()

    return number_reads, outfiles


def main(sample_id, fastq_files):

    logger.info("STARTING renamePE_samtoolsFASTQ.py")

    file_objects = []

    for fastq in fastq_files:

        logger.info("Processing file {}".format(fastq))

        logger.info("[{}] Guessing file compression".format(fastq))
        ftype = guess_file_compression(fastq)

        # This can guess the compression of gz, bz2 and zip. If it cannot
        # find the compression type, it tries to open a regular file.
        if ftype:
            logger.info("[{}] Found file compression: {}".format(fastq, ftype))
            file_objects.append(COPEN[ftype](fastq, "rt"))

        else:
            logger.info("[{}] File compression not found. Assuming an uncompressed file".format(fastq))
            file_objects.append(fastq)

    logger.info('Renaming fastq headers')
    number_reads, outfiles = formartFastqHeaders(sample_id, file_objects[0], file_objects[1])

    logger.info('{} read pairs were written in {} and {}. Compressing...'.format(number_reads, outfiles[0], outfiles[1]))

    # compress outfiles
    for file in outfiles:
        with open(file, 'rb') as f_in:
            f_out = gzip.open(file + '.gz', 'wb')
            f_out.writelines(f_in)
            f_out.close()
    logger.info('DONE')

    os.remove(outfiles[0])
    os.remove((outfiles[1]))


if __name__ == "__main__":
    main(SAMPLE_ID, FASTQ_PAIR)
