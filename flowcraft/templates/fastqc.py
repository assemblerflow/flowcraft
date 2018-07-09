#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to run FastQC on paired-end FastQ files.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``fastq_pair`` : *Pair of FastQ file paths*
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'``

Generated output
----------------

The generated output are output files that contain an object, usually a string.

- ``pair_{1,2}_data`` : File containing FastQC report at the nucleotide level\
    for each pair
    - e.g.: ``'pair_1_data'`` and ``'pair_2_data'``
- ``pair_{1,2}_summary``: File containing FastQC report for each category and\
    for each pair
    - e.g.: ``'pair_1_summary'`` and ``'pair_2_summary'``

Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "28032018"
__template__ = "fastqc-nf"

import os
import subprocess

from subprocess import PIPE
from os.path import exists, join

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


def __get_version_fastqc():

    try:

        cli = ["fastqc", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        version = stdout.strip().split()[1][1:].decode("utf8")

    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return {
        "program": "FastQC",
        "version": version,
    }


if __file__.endswith(".command.sh"):
    FASTQ_PAIR = '$fastq_pair'.split()
    ADAPTER_FILE = '$ad'
    CPUS = '$task.cpus'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))
    logger.debug("ADAPTER_FILE: {}".format(ADAPTER_FILE))
    logger.debug("CPUS: {}".format(CPUS))


def convert_adatpers(adapter_fasta):
    """Generates an adapter file for FastQC from a fasta file.

    The provided adapters file is assumed to be a simple fasta file with the
    adapter's name as header and the corresponding sequence::

        >TruSeq_Universal_Adapter
        AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
        >TruSeq_Adapter_Index 1
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG

    Parameters
    ----------
    adapter_fasta : str
        Path to Fasta file with adapter sequences.

    Returns
    -------
    adapter_out : str or None
        The path to the reformatted adapter file. Returns ``None`` if the
        adapters file does not exist or the path is incorrect.
    """

    adapter_out = "fastqc_adapters.tab"
    logger.debug("Setting output adapters file to: {}".format(adapter_out))

    try:

        with open(adapter_fasta) as fh, \
                open(adapter_out, "w") as adap_fh:

            for line in fh:
                if line.startswith(">"):

                    head = line[1:].strip()
                    # Get the next line with the sequence string
                    sequence = next(fh).strip()

                    adap_fh.write("{}\\t{}\\n".format(head, sequence))

        logger.info("Converted adapters file")

        return adapter_out

    # If an invalid adapters file is provided, return None.
    except FileNotFoundError:
        logger.warning("Could not find the provided adapters file: {}".format(
            adapter_fasta))
        return


@MainWrapper
def main(fastq_pair, adapter_file, cpus):
    """ Main executor of the fastq template.

    Parameters
    ----------
    fastq_pair : list
        Two element list containing the paired FastQ files.
    adapter_file : str
        Path to adapters file.
    cpus : int or str
        Number of cpu's that will be by FastQC.

    """

    logger.info("Starting fastqc")

    # If an adapter file was provided, convert it to FastQC format
    if os.path.exists(adapter_file):
        logger.info("Adapters file provided: {}".format(adapter_file))
        adapters = convert_adatpers(adapter_file)
    else:
        logger.info("Adapters file '{}' not provided or does not "
                    "exist".format(adapter_file))
        adapters = None

    # Setting command line for FastQC
    cli = [
        "fastqc",
        "--extract",
        "--nogroup",
        "--format",
        "fastq",
        "--threads",
        str(cpus)
    ]

    # Add adapters file to command line, if it exists
    if adapters:
        cli += ["--adapters", "{}".format(adapters)]

    # Add FastQ files at the end of command line
    cli += fastq_pair

    logger.debug("Running fastqc subprocess with command: {}".format(cli))

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE, shell=False)
    stdout, stderr = p.communicate()

    # Attempt to decode STDERR output from bytes. If unsuccessful, coerce to
    # string
    try:
        stderr = stderr.decode("utf8")
    except (UnicodeDecodeError, AttributeError):
        stderr = str(stderr)

    logger.info("Finished fastqc subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished fastqc subprocesswith STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished fastqc with return code: {}".format(
        p.returncode))

    logger.info("Checking if FastQC output was correctly generated")
    # Check if the FastQC output was correctly generated.
    with open(".status", "w") as status_fh:
        for fastq in fastq_pair:
            fpath = join(fastq.rsplit(".", 2)[0] + "_fastqc",
                         "fastqc_data.txt")
            logger.debug("Checking path: {}".format(fpath))
            # If the FastQC output does not exist, pass the STDERR to
            # the output status channel and exit
            if not exists(fpath):
                logger.warning("Path does not exist: {}".format(fpath))
                status_fh.write("fail")
                return

            logger.debug("Found path: {}".format(fpath))

        # If the output directories exist, write 'pass' to the output status
        # channel
            status_fh.write("pass")

    logger.info("Retrieving relevant FastQC output files")

    # Both FastQC have been correctly executed. Get the relevant FastQC
    # output files for the output channel
    for i, fastq in enumerate(fastq_pair):
        # Get results for each pair
        fastqc_dir = fastq.rsplit(".", 2)[0] + "_fastqc"

        summary_file = join(fastqc_dir, "summary.txt")
        logger.debug("Retrieving summary file: {}".format(summary_file))
        fastqc_data_file = join(fastqc_dir, "fastqc_data.txt")
        logger.debug("Retrieving data file: {}".format(fastqc_data_file))

        # Rename output files to a file name that is easier to handle in the
        # output channel
        os.rename(fastqc_data_file, "pair_{}_data".format(i + 1))
        os.rename(summary_file, "pair_{}_summary".format(i + 1))


if __name__ == "__main__":

    main(FASTQ_PAIR, ADAPTER_FILE, CPUS)
