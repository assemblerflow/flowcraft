#!/usr/bin/env python3


"""
Purpose
-------

This module is intended to generate a json output for mapping results that
can be imported in pATLAS.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``depth_file`` : String with the name of the mash screen output file.
    - e.g.: ``'samtoolsDepthOutput_sampleA.txt'``
- ``json_dict`` : the file that contains the dictionary with keys and values for
        accessions and their respective lengths.
    - e.g.: ``'reads_sample_result_length.json'``
- ``cutoff`` : The cutoff used to trim the unwanted matches for the minimum
        coverage results from mapping. This value may range between 0 and 1.
    - e.g.: ``0.6``


Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "20022018"
__template__ = "mapping2json-nf"

import os
import json

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    DEPTH_TXT = '$depthFile'
    JSON_LENGTH = '$lengthJson'
    CUTOFF = '$params.cov_cutoff'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("DEPTH_TXT: {}".format(DEPTH_TXT))
    logger.debug("JSON_LENGHT: {}".format(JSON_LENGTH))
    logger.debug("CUTOFF: {}".format(CUTOFF))


def depthfilereader(depth_file, plasmid_length, cutoff):
    '''
    Function that parse samtools depth file and creates 3 dictionaries that
    will be useful to make the outputs of this script, both the tabular file
    and the json file that may be imported by pATLAS

    Parameters
    ----------
    depth_file: str
        the path to depth file for each sample
    plasmid_length: dict
        a dictionary that stores length of all plasmids in fasta given as input
    cutoff: str
        the cutoff used to trim the unwanted matches for the minimum coverage
        results from mapping. This is then converted into a float within this
        function in order to compare with the value returned from the perc_value_per_ref.

    Returns
    -------
    percentage_basescovered: dict
            stores the percentage of the total sequence of a
            reference/accession (plasmid) in a dictionary
    '''
    depth_dic_coverage = {}
    for line in depth_file:
        tab_split = line.split()    # split by any white space
        reference = "_".join(tab_split[0].strip().split("_")[0:3])  # store
        # only the gi for the reference
        position = tab_split[1]
        numreadsalign = float(tab_split[2].rstrip())
        if reference not in depth_dic_coverage:
            depth_dic_coverage[reference] = {}
        depth_dic_coverage[reference][position] = numreadsalign

    percentage_basescovered = {}
    for ref in depth_dic_coverage:
        # calculates the percentage value per each reference
        perc_value_per_ref = float(len(depth_dic_coverage[ref])) / \
                                       float(plasmid_length[ref])
        # checks if percentage value is higher or equal to the cutoff defined
        if perc_value_per_ref >= float(cutoff):
            percentage_basescovered[ref] = perc_value_per_ref

    return percentage_basescovered

@MainWrapper
def main(depth_file, json_dict, cutoff):
    '''
    Function that handles the inputs required to parse depth files from bowtie
    and dumps a dict to a json file that can be imported into pATLAS.

    Parameters
    ----------
    depth_file: str
         the path to depth file for each sample
    json_dict: str
        the file that contains the dictionary with keys and values for accessions
        and their respective lengths
    cutoff: str
        the cutoff used to trim the unwanted matches for the minimum coverage
        results from mapping. This value may range between 0 and 1.


    '''

    # check for the appropriate value for the cutoff value for coverage results
    try:
        cutoff_val = float(cutoff)
    except ValueError:
        logger.error("Cutoff value should be a string such as: '0.6'. "
                     "The outputted value: {}. Make sure to provide an "
                     "appropriate value for --cov_cutoff".format(cutoff))

    # loads dict from file, this file is provided in docker image

    plasmid_length = json.load(open(json_dict))

    # read depth file
    depth_file_reader = open(depth_file)

    # first reads the depth file and generates dictionaries to handle the input
    # to a simpler format
    logger.info("Reading depth file and creating dictionary to dump")
    percentage_basescovered = depthfilereader(depth_file_reader, plasmid_length,
                                              cutoff_val)

    # then dump do file
    output_json = open("{}_mapping.json".format(depth_file), "w")
    logger.info("Dumping to {}".format("{}_mapping.json".format(depth_file)))
    output_json.write(json.dumps(percentage_basescovered))
    output_json.close()


if __name__ == "__main__":

    main(DEPTH_TXT, JSON_LENGTH, CUTOFF)
