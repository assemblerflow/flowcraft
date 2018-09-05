#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to generate a json output for mash screen results that
can be imported in pATLAS.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``mash_output`` : String with the name of the mash screen output file.
    - e.g.: ``'sortedMashScreenResults_SampleA.txt'``


Code documentation
------------------

"""

__version__ = "1.1.0"
__build__ = "04072018"
__template__ = "mashscreen2json-nf"

from statistics import median
import os
import json

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    MASH_TXT = '$mashtxt'
    SAMPLE_ID = '$sample_id'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("MASH_TXT: {}".format(MASH_TXT))
    logger.debug("SAMPLE_ID: {}".format(MASH_TXT))

@MainWrapper
def main(mash_output, sample_id):
    '''
    converts top results from mash screen txt output to json format

    Parameters
    ----------
    mash_output: str
        this is a string that stores the path to this file, i.e, the name of
        the file
    sample_id: str
        sample name

    '''
    logger.info("Reading file : {}".format(mash_output))
    read_mash_output = open(mash_output)

    dic = {}
    median_list = []
    filtered_dic = {}

    logger.info("Generating dictionary and list to pre-process the final json")
    for line in read_mash_output:
        tab_split = line.split("\t")
        identity = tab_split[0]
        # shared_hashes = tab_split[1]
        median_multiplicity = tab_split[2]
        # p_value = tab_split[3]
        query_id = tab_split[4]
        # query-comment should not exist here and it is irrelevant

        # here identity is what in fact interests to report to json but
        # median_multiplicity also is important since it gives an rough
        # estimation of the coverage depth for each plasmid.
        # Plasmids should have higher coverage depth due to their increased
        # copy number in relation to the chromosome.
        dic[query_id] = [identity, median_multiplicity]
        median_list.append(float(median_multiplicity))

    output_json = open(" ".join(mash_output.split(".")[:-1]) + ".json", "w")

    # median cutoff is twice the median of all median_multiplicity values
    # reported by mash screen. In the case of plasmids, since the database
    # has 9k entries and reads shouldn't have that many sequences it seems ok...
    if len(median_list) > 0:
        # this statement assures that median_list has indeed any entries
        median_cutoff = median(median_list)
        logger.info("Generating final json to dump to a file")
        for k, v in dic.items():
            # estimated copy number
            copy_number = int(float(v[1]) / median_cutoff)
            # assure that plasmid as at least twice the median coverage depth
            if float(v[1]) > median_cutoff:
                filtered_dic["_".join(k.split("_")[0:3])] = [
                    round(float(v[0]),2),
                    copy_number
                ]
        logger.info(
            "Exported dictionary has {} entries".format(len(filtered_dic)))
    else:
        # if no entries were found raise an error
        logger.error("No matches were found using mash screen for the queried reads")

    output_json.write(json.dumps(filtered_dic))
    output_json.close()

    json_dic = {
        "tableRow": [{
            "sample": sample_id,
            "data": [{
                "header": "Mash Screen",
                "table": "plasmids",
                "patlas_mashscreen": filtered_dic,
                "value": len(filtered_dic)
            }]
        }],
    }

    with open(".report.json", "w") as json_report:
        json_report.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == "__main__":

    main(MASH_TXT, SAMPLE_ID)
