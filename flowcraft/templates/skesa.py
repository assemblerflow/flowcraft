#!/usr/bin/env python3

"""
Purpose
-------

This module is intended execute Skesa on paired-end FastQ files.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fastq_pair`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'``

Generated output
----------------

- ``${sample_id}_*.assembly.fasta`` : Main output of skesawith the assembly
    - e.g.: ``sample_1_skesa.fasta``

Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "16012018"
__template__ = "skesa-nf"

import os
import re
import subprocess

from subprocess import PIPE

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


def __get_version_skesa():

    try:

        cli = ["skesa", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        _, err = p.communicate()

        try:
            version = re.search("v((\\..*))-", err.decode("utf8")).group(1)
        except AttributeError:
            version = "undefined"

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
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))


@MainWrapper
def main(sample_id, fastq_pair):
    """Main executor of the skesa template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    fastq_pair : list
        Two element list containing the paired FastQ files.
    """

    logger.info("Starting skesa")

    # Determine output file
    if "_trim." in fastq_pair[0]:
        sample_id += "_trim"
    version = __get_version_skesa()["version"]
    output_file = "{}_skesa{}.fasta".format(sample_id, version.replace(".", ""))

    cli = [
        "skesa",
        "--fastq",
        "{},{}".format(fastq_pair[0], fastq_pair[1]),
        "--gz",
        "--use_paired_ends",
        "--cores",
        "${task.cpus}"
    ]

    logger.debug("Running Skesa subprocess with command: {}".format(cli))

    with open(output_file, "w") as fh:
        p = subprocess.Popen(cli, stdout=fh, stderr=PIPE)
    stdout, stderr = p.communicate()

    # Attempt to decode STDERR output from bytes. If unsuccessful, coerce to
    # string
    try:
        stderr = stderr.decode("utf8")
        stdout = stdout.decode("utf8")
    except (UnicodeDecodeError, AttributeError):
        stderr = str(stderr)
        stdout = str(stdout)

    logger.info("Finished Skesa subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished Skesa subprocesswith STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished Skesa with return code: {}".format(
        p.returncode))

    with open(".status", "w") as fh:
        if p.returncode != 0:
            fh.write("error")
            raise SystemExit(p.returncode)
        else:
            fh.write("pass")


if __name__ == '__main__':

    main(SAMPLE_ID, FASTQ_PAIR)
