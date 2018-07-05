#!/usr/bin/env python3

"""
Purpose
-------

This module is intended execute metaSpades on paired-end FastQ files.

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

__version__ = "1.0.1"
__build__ = "16012018"
__template__ = "metaspades-nf"

import os
import subprocess

from subprocess import PIPE

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


def __get_version_spades():

    try:

        cli = ["metaspades.py", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        version = stdout.strip().split()[-1][1:].decode("utf8")

    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return {
        "program": "metaSPAdes",
        "version": version,
    }


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    FASTQ_PAIR = '$fastq_pair'.split()
    MAX_LEN = int('$max_len'.strip())
    KMERS = '$kmers'.strip()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))
    logger.debug("MAX_LEN: {}".format(MAX_LEN))
    logger.debug("KMERS: {}".format(KMERS))


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


@MainWrapper
def main(sample_id, fastq_pair, max_len, kmer):
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

    """

    logger.info("Starting spades")

    logger.info("Setting SPAdes kmers")
    kmers = set_kmers(kmer, max_len)
    logger.info("SPAdes kmers set to: {}".format(kmers))

    cli = [
        "metaspades.py",
        "--only-assembler",
        "--threads",
        "$task.cpus",
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

    logger.debug("Running metaSPAdes subprocess with command: {}".format(cli))

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

    logger.info("Finished metaSPAdes subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished metaSPAdes subprocesswith STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished metaSPAdes with return code: {}".format(
        p.returncode))

    with open(".status", "w") as fh:
        if p.returncode != 0:
            fh.write("error")
            return
        else:
            fh.write("pass")

    # Change the default contigs.fasta assembly name to a more informative one
    if "_trim." in fastq_pair[0]:
        sample_id += "_trim"

    assembly_file = "{}_metaspades.fasta".format(
        sample_id)
    os.rename("contigs.fasta", assembly_file)
    logger.info("Setting main assembly file to: {}".format(assembly_file))


if __name__ == '__main__':

    main(SAMPLE_ID, FASTQ_PAIR, MAX_LEN, KMERS)
