#!/usr/bin/env python3

"""
Purpose
-------

This module is intended execute Spades on paired-end FastQ files.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fastq_pair`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'``
- ``kmers`` : Setting for Spades kmers. Can be either ``'auto'``, \
    ``'default'`` or a user provided list.
    - e.g.: ``'auto'`` or ``'default'`` or ``'55 77 99 113 127'``
- ``opts`` : List of options for spades execution.
    1. The minimum number of reads to consider an edge in the de Bruijn \
    graph during the assembly.
        - e.g.: ``'5'``
    2. Minimum contigs k-mer coverage.
        - e.g.: ``['2' '2']``
- ``clear`` : If 'true', remove the input fastq files at the end of the
    component run, IF THE FILES ARE IN THE WORK DIRECTORY

Generated output
----------------

- ``contigs.fasta`` : Main output of spades with the assembly
    - e.g.: ``contigs.fasta``
- ``spades_status`` :  Stores the status of the spades run. If it was \
    successfully executed, it stores ``'pass'``. Otherwise, it stores the\
    ``STDERR`` message.
    - e.g.: ``'pass'``

Code documentation
------------------

"""

__version__ = "1.0.2"
__build__ = "29062018"
__template__ = "spades-nf"

import os
import sys
import re
import subprocess

from subprocess import PIPE

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


def __get_version_spades():

    try:

        cli = ["spades.py", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        version = stdout.strip().split()[-1][1:].decode("utf8")

    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return {
        "program": "SPAdes",
        "version": version,
    }


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    FASTQ_PAIR = '$fastq_pair'.split()
    MAX_LEN = int('$max_len'.strip())
    KMERS = '$kmers'.strip()
    CLEAR = '$clear'
    OPTS = [x.strip() for x in '$opts'.strip("[]").split(",")]
    CLEAR = '$clear'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))
    logger.debug("MAX_LEN: {}".format(MAX_LEN))
    logger.debug("KMERS: {}".format(KMERS))
    logger.debug("OPTS: {}".format(OPTS))
    logger.debug("CLEAR: {}".format(CLEAR))


def set_kmers(kmer_opt, max_read_len):
    """Returns a kmer list based on the provided kmer option and max read len.

    Parameters
    ----------
    kmer_opt : str
        The k-mer option. Can be either ``'auto'``, ``'default'`` or a
        sequence of space separated integers, ``'23, 45, 67'``.
    max_read_len : int
        The maximum read length of the current sample.

    Returns
    -------
    kmers : list
        List of k-mer values that will be provided to Spades.

    """

    logger.debug("Kmer option set to: {}".format(kmer_opt))

    # Check if kmer option is set to auto
    if kmer_opt == "auto":

        if max_read_len >= 175:
            kmers = [55, 77, 99, 113, 127]
        else:
            kmers = [21, 33, 55, 67, 77]

        logger.debug("Kmer range automatically selected based on max read"
                     "length of {}: {}".format(max_read_len, kmers))

    # Check if manual kmers were specified
    elif len(kmer_opt.split()) > 1:

        kmers = kmer_opt.split()
        logger.debug("Kmer range manually set to: {}".format(kmers))

    else:

        kmers = []
        logger.debug("Kmer range set to empty (will be automatically "
                     "determined by SPAdes")

    return kmers


def clean_up(fastq):
    """
    Cleans the temporary fastq files. If they are symlinks, the link
    source is removed

    Parameters
    ----------
    fastq : list
        List of fastq files.
    """

    for fq in fastq:
        # Get real path of fastq files, following symlinks
        rp = os.path.realpath(fq)
        logger.debug("Removing temporary fastq file path: {}".format(rp))
        if re.match(".*/work/.{2}/.{30}/.*", rp):
            os.remove(rp)


@MainWrapper
def main(sample_id, fastq_pair, max_len, kmer, opts, clear):
    """Main executor of the spades template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    fastq_pair : list
        Two element list containing the paired FastQ files.
    max_len : int
        Maximum read length. This value is determined in
        :py:class:`templates.integrity_coverage`
    kmer : str
        Can be either ``'auto'``, ``'default'`` or a
        sequence of space separated integers, ``'23, 45, 67'``.
    opts : List of options for spades execution. See above.
    clear : str
        Can be either 'true' or 'false'. If 'true', the input fastq files will
        be removed at the end of the run, IF they are in the working directory

    """

    logger.info("Starting spades")

    min_coverage, min_kmer_coverage = opts

    logger.info("Setting SPAdes kmers")
    kmers = set_kmers(kmer, max_len)
    logger.info("SPAdes kmers set to: {}".format(kmers))

    cli = [
        "spades.py",
        "--careful",
        "--only-assembler",
        "--threads",
        "$task.cpus",
        "--cov-cutoff",
        min_coverage,
        "-o",
        "."
    ]

    # Add kmers, if any were specified
    if kmers:
        cli += ["-k {}".format(",".join([str(x) for x in kmers]))]

    # Add FastQ files
    cli += [
        "-1",
        fastq_pair[0],
        "-2",
        fastq_pair[1]
    ]

    logger.debug("Running SPAdes subprocess with command: {}".format(cli))

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # Attempt to decode STDERR output from bytes. If unsuccessful, coerce to
    # string
    try:
        stderr = stderr.decode("utf8")
        stdout = stdout.decode("utf8")
    except (UnicodeDecodeError, AttributeError):
        stderr = str(stderr)
        stdout = str(stdout)

    logger.info("Finished SPAdes subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished SPAdes subprocess with STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished SPAdes with return code: {}".format(
        p.returncode))

    with open(".status", "w") as fh:
        if p.returncode != 0:
            fh.write("error")
            sys.exit(p.returncode)
        else:
            fh.write("pass")

    # Change the default contigs.fasta assembly name to a more informative one
    if "_trim." in fastq_pair[0]:
        sample_id += "_trim"
    # Get spades version for output name
    info = __get_version_spades()

    assembly_file = "{}_spades{}.fasta".format(
        sample_id, info["version"].replace(".", ""))
    os.rename("contigs.fasta", assembly_file)
    logger.info("Setting main assembly file to: {}".format(assembly_file))

    # Remove input fastq files when clear option is specified.
    # Only remove temporary input when the expected output exists.
    if clear == "true" and os.path.exists("contigs.fasta"):
        clean_up(fastq_pair)


if __name__ == '__main__':

    main(SAMPLE_ID, FASTQ_PAIR, MAX_LEN, KMERS, OPTS, CLEAR)
