#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to generate a json output from the consensus results from
all the approaches available through options (mapping, assembly, mash screen)

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``mapping_json`` : String with the name of the json file with mapping results.
    - e.g.: ``'mapping_SampleA.json'``
- ``dist_json`` : String with the name of the json file with mash dist results.
    - e.g.: ``'mash_dist_SampleA.json'``
- ``screen_json`` : String with the name of the json file with mash screen results.
    - e.g.: ``'mash_screen_sampleA.json'``


Code documentation
------------------

"""

__version__ = "0.1.0"
__build__ = "24022018"
__template__ = "pATLAS_consensus_json-nf"

import os
import json

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    LIST_OF_FILES = '$infile_list'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("LIST_OF_FILES: {}".format(LIST_OF_FILES))


@MainWrapper
def main(list_of_jsons):
    """

    Parameters
    ----------
    list_of_jsons: list
        A list of files provided by fullConsensus process provided by nextflow

    """

    # first lets gather a collection of the input and their corresponding dicts
    file_correspondence = {}

    for infile in list_of_jsons:
        file_dict = json.load(open(infile))
        file_correspondence[infile] = file_dict

    json_dict = {}
    for accession in list(file_correspondence.values())[0]:
        if all([True if accession in f_dict else False
                for f_dict in file_correspondence.values()]):
            accession_dict = {}
            for infile in file_correspondence.keys():
                accession_dict[infile] = file_correspondence[infile][accession]

            json_dict[accession] = accession_dict

    out_file = open("consensus_{}.json".format(
        list_of_jsons[0].split(".")[0].split("_")[-1]), "w")

    out_file.write(json.dumps(json_dict))
    out_file.close()

    json_dic = {
        "patlas_mashscreen": json_dict
        # TODO add information for report webapp
    }

    with open(".report.json", "w") as json_report:
        json_report.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == "__main__":
    main(LIST_OF_FILES)