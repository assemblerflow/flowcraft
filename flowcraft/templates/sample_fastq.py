#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to sub-sample FastQ files to a certain coverage, based
on the expected genome size.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fastq_pair`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'``
- ``gsize`` : *Expected genome size*
    - e.g.: ``'2.5'``
- ``depth`` : Maximum depth threshold above which the subsampling will be
    performed.
    - e.g.: ``100``
- ``clear`` : If 'true', remove the input fastq files at the end of the
    component run, IF THE FILES ARE IN THE WORK DIRECTORY

Generated output
----------------

- ``*_ss.fq.gz`` : Subsample fastq reads
    - e.g.: ``sampleA_ss.fq.gz``

Code documentation
------------------

"""

__version__ = "1.0.0"
__build__ = "30072018"
__template__ = "sample_fastq-nf"

import os
import re
import json
import subprocess

from os.path import basename

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    FASTQ_PAIR = '$fastq_pair'.split()
    GSIZE = float('$gsize'.strip())
    DEPTH = float('$depth'.strip())
    CLEAR = '$clear'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))
    logger.debug("GENOME_SIZE: {}".format(GSIZE))
    logger.debug("DEPTH: {}".format(DEPTH))
    logger.debug("CLEAR: {}".format(CLEAR))


def __get_version_spades():

    try:

        cli = ["seqtk"]
        p = subprocess.Popen(cli, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        _, stderr = p.communicate()

        _version = stderr.splitlines()[2]
        try:
            version = re.match(
                "Version: (.*)", _version.decode("utf8")).group(1)
        except AttributeError:
            version = "undefined"

    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return {
        "program": "seqtk",
        "version": version,
    }


@MainWrapper
def main(sample_id, fastq_pair, genome_size, depth, clear):

    genome_size = genome_size
    target_depth = depth
    p1 = fastq_pair[0]
    p2 = fastq_pair[1]
    bn1 = ".".join(basename(p1).split('.')[:-2])
    bn2 = ".".join(basename(p2).split('.')[:-2])

    R1_fqchk = subprocess.Popen(['seqtk', 'fqchk', p1], stdout=subprocess.PIPE)
    R1_stdout, R1_stderr = R1_fqchk.communicate()
    B_P1 = int(R1_stdout.splitlines()[2].split()[1])
    logger.debug("Bases p1: {}".format(B_P1))

    R2_fqchk = subprocess.Popen(['seqtk', 'fqchk', p2], stdout=subprocess.PIPE)
    R2_stdout, R2_stderr = R2_fqchk.communicate()
    B_P2= int(R2_stdout.splitlines()[2].split()[1])
    logger.debug("Bases p2: {}".format(B_P2))

    estimated_coverage = (B_P1 + B_P2) / (genome_size * 1E6)
    logger.debug("Estimated coverage: {}".format(estimated_coverage))
    ratio = target_depth/estimated_coverage
    logger.debug("Estimated ration: {}".format(ratio))

    if ratio < 1:
        # print ("Writing R1.fq.gz")
        ps = subprocess.Popen(('seqtk', 'sample', '-s100', p1, str(ratio)),
                              stdout=subprocess.PIPE)
        with open('{}_ss.fq.gz'.format(bn1), 'w') as outfile:
            subprocess.Popen(('gzip', '--fast', '-c'),
                             stdin=ps.stdout, stdout=outfile )
        ps.wait()

        # print ("Writing R2.fq.gz")
        ps = subprocess.Popen(('seqtk', 'sample', '-s100', p2, str(ratio)),
                              stdout=subprocess.PIPE)
        with open('{}_ss.fq.gz'.format(bn2), 'w') as outfile:
            subprocess.Popen(('gzip', '--fast', '-c'),
                             stdin=ps.stdout, stdout=outfile)
        ps.wait()

        if clear == "true":
            # Get real path of the symlink
            for fq in [p1, p2]:
                rp = os.path.realpath(fq)
                print("removing temporary fastq file path: {}".format(rp))
                # remove only when the file is in the work directory
                if re.match(".*/work/.{2}/.{30}/.*", rp):
                    os.remove(rp)

    else:
        os.symlink(p1, "{}._ss.fq.gz".format(bn1))
        os.symlink(p2, "{}._ss.fq.gz".format(bn2))

    # Record the original estimated coverage
    with open(".report.json", "w") as fh:
        json_dic = {
            "tableRow": [
                {
                    "sample": sample_id,
                    "data": [{
                        "header": "Coverage",
                        "value": round(estimated_coverage, 1),
                        "table": "qc"
                    }]
                 }
            ]
        }
        fh.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == "__main__":
        main(SAMPLE_ID, FASTQ_PAIR, GSIZE, DEPTH, CLEAR)
