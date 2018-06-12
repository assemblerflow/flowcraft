#!/usr/bin/env python3


"""
Purpose
-------

This module is intended to generate a json output for mash dist results that
can be imported in pATLAS.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``mash_output`` : String with the name of the mash screen output file.
    - e.g.: ``'fastaFileA_mashdist.txt'``


Code documentation
------------------

"""

__version__ = "1.2.0"
__build__ = "17052018"
__template__ = "mashsdist2json-nf"

import os
import json

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    MASH_TXT = '$mashtxt'
    HASH_CUTOFF = '$params.shared_hashes'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("MASH_TXT: {}".format(MASH_TXT))
    logger.debug("HASH_CUTOFF: {}".format(HASH_CUTOFF))


def send_to_output(master_dict, mash_output):
    """Send dictionary to output json file
    This function sends master_dict dictionary to a json file if master_dict is
    populated with entries, otherwise it won't create the file

    Parameters
    ----------
    master_dict: dict
        dictionary that stores all entries for a specific query sequence
        in multi-fasta given to mash dist as input against patlas database
    last_seq: str
        string that stores the last sequence that was parsed before writing to
        file and therefore after the change of query sequence between different
        rows on the input file
    mash_output: str
        the name/path of input file to main function, i.e., the name/path of
        the mash dist output txt file.

    Returns
    -------

    """
    # create a new file only if master_dict is populated
    if master_dict:
        out_file = open("{}.json".format(
            "".join(mash_output.split(".")[0])), "w")
        out_file.write(json.dumps(master_dict))
        out_file.close()


@MainWrapper
def main(mash_output, hash_cutoff):
    '''
    Main function that allows to dump a mash dist txt file to a json file

    Parameters
    ----------
    mash_output: str
        A string with the input file.

    '''
    # out_file = open(".".join(mash_output.split(".")[:-1]) + ".json", "w")
    input_f = open(mash_output, "r")

    master_dict = {}
    # used to store the last sequence to be parsed (useful for multifasta)
    last_seq = ""
    counter = 0

    for line in input_f:

        tab_split = line.split("\t")
        current_seq = tab_split[1].strip()
        ref_accession = "_".join(tab_split[0].strip().split("_")[0:3])
        mash_dist = tab_split[2].strip()
        hashes_list = tab_split[-1].strip().split("/")

        # creates a percentage of the shared hashes between the sample and the
        # reference
        perc_hashes = float(hashes_list[0]) / float(hashes_list[1])

        # assures that only the hashes with a given shared percentage are
        # reported to json file
        if perc_hashes > float(hash_cutoff):

            master_dict[ref_accession] = [1 - float(mash_dist), perc_hashes,
                                              current_seq]

    # assures that file is closed in last iteration of the loop
    send_to_output(master_dict, mash_output)


if __name__ == "__main__":

    main(MASH_TXT, HASH_CUTOFF)
