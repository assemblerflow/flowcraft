#!/usr/bin/env python3

"""
Purpose
-------

This module intends to split a multifasta file into seperate fasta files.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fasta`` : A fasta file path.
    - e.g.: ``'SampleA.fasta'``
- ``min_contig_size`` : A minimum contig length
    - e.g.: ``'1000'``

Generated output
----------------

-  A fasta file per contig (given the minimum contig size
"""

__version__ = "0.0.2"
__build__ = "08082018"
__template__ = "split_assembly-nf"

import os
from itertools import groupby

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLY = '$assembly'
    MIN_SIZE = int('$min_contig_size'.strip())
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("MIN_SIZE: {}".format(MIN_SIZE))



@MainWrapper
def main(sample_id, assembly, min_size):
    """Main executor of the split_fasta template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    assembly : list
        Assembly file.
    min_size : int
        Minimum contig size.

    """

    logger.info("Starting script")

    f_open = open(assembly, "rU")

    entry = (x[1] for x in groupby(f_open, lambda line: line[0] == ">"))

    success = 0

    for header in entry:
        headerStr = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in entry.__next__())
        if len(seq) >= min_size:
            with open(sample_id + '_' + headerStr.replace(" ","_").replace("=","_") + '.fasta', "w") as output_file:
                output_file.write(">" + sample_id + "_" + headerStr.replace(" ","_").replace("=","_") + "\\n" + seq + "\\n")
                success += 1

    f_open.close()

    logger.info("{} sequences sucessfully splitted.".format(success))


if __name__ == '__main__':

    main(SAMPLE_ID, ASSEMBLY, MIN_SIZE)