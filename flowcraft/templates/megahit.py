#!/usr/bin/env python3

"""
Purpose
-------

This module is intended execute megahit on paired-end FastQ files.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fastq_pair`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'``
- ``kmers`` : Setting for megahit kmers. Can be either ``'auto'``, \
    ``'default'`` or a user provided list. All must be odd, in the range 15-255, increment <= 28
    - e.g.: ``'auto'`` or ``'default'`` or ``'55 77 99 113 127'``

Generated output
----------------

- ``contigs.fa`` : Main output of megahit with the assembly
    - e.g.: ``contigs.fa``
- ``megahit_status`` :  Stores the status of the megahit run. If it was \
    successfully executed, it stores ``'pass'``. Otherwise, it stores the\
    ``STDERR`` message.
    - e.g.: ``'pass'``

Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "26042018"
__template__ = "megahit-nf"

import os
import subprocess

from subprocess import PIPE

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


def is_odd(k_mer):
    for i in k_mer:
        if i % 2 != 0:
            return True
    return False


def __get_version_megahit():

    try:

        cli = ["megahit", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        version = stdout.strip().split()[-1][1:].decode("utf8")

    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return {
        "program": "megahit",
        "version": version,
    }


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    FASTQ_PAIR = '$fastq_pair'.split()
    MAX_LEN = int('$max_len'.strip())
    KMERS = '$kmers'.strip()
    MEM = '$task.memory'
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
        List of k-mer values that will be provided to megahit.

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

    # Check if manual k-mers were specified
    elif len(kmer_opt.split()) > 1:

        kmers = kmer_opt.split()
        if kmers[0]<15 or kmers[-1]>255 or is_odd(kmers):
            kmers = []
            logger.debug("Kmer out of range or with even numbers"
                         "(will be automatically determined by megahit")
        else:
            logger.debug("Kmer range manually set to: {}".format(kmers))

    else:

        kmers = []
        logger.debug("Kmer range set to empty (will be automatically "
                     "determined by megahit")

    return kmers


def fix_contig_names(asseembly_path):
    """Removes whitespace from the assembly contig names

    Parameters
    ----------
    asseembly_path : path to assembly file

    Returns
    -------
    str:
        Path to new assembly file with fixed contig names
    """

    fixed_assembly = "fixed_assembly.fa"

    with open(asseembly_path) as in_hf, open(fixed_assembly, "w") as ou_fh:

        for line in in_hf:

            if line.startswith(">"):
                fixed_line = line.replace(" ", "_")
                ou_fh.write(fixed_line)
            else:
                ou_fh.write(line)

    return fixed_assembly


@MainWrapper
def main(sample_id, fastq_pair, max_len, kmer, mem):
    """Main executor of the megahit template.

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

    logger.info("Starting megahit")

    logger.info("Setting megahit kmers")
    kmers = set_kmers(kmer, max_len)
    logger.info("megahit kmers set to: {}".format(kmers))

    mem_bytes = int(mem.replace(" GB", "")) * 1073741824

    cli = [
        "megahit",
        "--num-cpu-threads",
        "$task.cpus",
        "--memory",
        str(mem_bytes),
        "-o",
        "megahit"
    ]

    # Add kmers, if any were specified
    if kmers:
        cli += [
            "--k-list",
            ",".join([str(x) for x in kmers])
        ]

    # Add FastQ files
    cli += [
        "-1",
        fastq_pair[0],
        "-2",
        fastq_pair[1]
    ]

    logger.debug("Running megahit subprocess with command: {}".format(cli))

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

    logger.info("Finished megahit subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished megahit subprocesswith STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished megahit with return code: {}".format(
        p.returncode))

    with open(".status", "w") as fh:
        if p.returncode != 0:
            fh.write("error")
            return
        else:
            fh.write("pass")

    assembly_path = "megahit/final.contigs.fa"
    fixed_assembly = fix_contig_names(assembly_path)

    # Change the default final.contigs.fa assembly name to a more informative
    #  one
    if "_trim." in fastq_pair[0]:
        sample_id += "_trim"
    # Get megahit version for output name
    info = __get_version_megahit()

    assembly_file = "{}_megahit{}.fasta".format(
        sample_id, info["version"].replace(".", ""))
    os.rename(fixed_assembly, assembly_file)
    logger.info("Setting main assembly file to: {}".format(assembly_file))


if __name__ == '__main__':

    main(SAMPLE_ID, FASTQ_PAIR, MAX_LEN, KMERS, MEM)
