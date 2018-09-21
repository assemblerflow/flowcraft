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

__version__ = "1.1.0"
__build__ = "04072018"
__template__ = "mapping2json-nf"

import os
import json
import sys
from pympler.asizeof import asizeof

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    DEPTH_TXT = '$depthFile'
    JSON_LENGTH = '$lengthJson'
    CUTOFF = '$cov_cutoff'
    SAMPLE_ID = '$sample_id'
else:
    DEPTH_TXT = sys.argv[1]
    JSON_LENGTH = sys.argv[2]
    CUTOFF = sys.argv[3]
    SAMPLE_ID = sys.argv[4]

logger.debug("List of arguments given: {}".format([
    DEPTH_TXT,
    JSON_LENGTH,
    CUTOFF,
    SAMPLE_ID
]))

# check if all variables are assigned
if DEPTH_TXT and JSON_LENGTH and SAMPLE_ID and CUTOFF:
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("DEPTH_TXT: {}".format(DEPTH_TXT))
    logger.debug("JSON_LENGHT: {}".format(JSON_LENGTH))
    logger.debug("CUTOFF: {}".format(CUTOFF))
else:
    logger.error("Args should be given to this template, either from sys.argv"
                 " or through nextflow variables")


def depth_file_reader(depth_file):
    """
    Function that parse samtools depth file and creates 3 dictionaries that
    will be useful to make the outputs of this script, both the tabular file
    and the json file that may be imported by pATLAS

    Parameters
    ----------
    depth_file: textIO
        the path to depth file for each sample

    Returns
    -------
    depth_dic_coverage: dict
            dictionary with the coverage per position for each plasmid
    """

    # dict to store the mean coverage for each reference
    depth_dic_coverage = {}

    for line in depth_file:
        tab_split = line.split()  # split by any white space
        reference = "_".join(tab_split[0].strip().split("_")[0:3])  # store
        # only the gi for the reference
        position = tab_split[1]
        num_reads_align = float(tab_split[2].rstrip())

        if reference not in depth_dic_coverage:
            depth_dic_coverage[reference] = {}

        depth_dic_coverage[reference][position] = num_reads_align

    logger.info("Finished parsing depth file.")
    depth_file.close()

    logger.debug("Size of dict_cov: {} kb".format(
        asizeof(depth_dic_coverage)/1024))

    return depth_dic_coverage


def generate_jsons(depth_dic_coverage, plasmid_length, cutoff):
    """

    Parameters
    ----------
    depth_dic_coverage: dict
         dictionary with the coverage per position for each plasmid

    Returns
    -------
    percentage_bases_covered: dict
    dict_cov:  dict

    """

    # initializes the dictionary with the mean coverage results per plasmid
    percentage_bases_covered = {}
    # dict to store coverage results for a given interval of points
    dict_cov = {}

    for ref in depth_dic_coverage:
        # calculates the percentage value per each reference
        perc_value_per_ref = float(len(depth_dic_coverage[ref])) / \
            float(plasmid_length[ref])
        # checks if percentage value is higher or equal to the cutoff defined
        if perc_value_per_ref >= cutoff:
            percentage_bases_covered[ref] = round(perc_value_per_ref, 2)

            # starts parser to get the array with the coverage for all the
            # positions
            # first, sets the interval for the reference being parsed
            interval = round(int(plasmid_length[ref]) * 0.01,
                             ndigits=0)

            # if the sequence is smaller than 100 bp, which shouldn't happen
            # anyway
            if interval < 1:
                interval = 1

            # starts dict cov for the reference
            dict_cov[ref] = {
                "length": int(plasmid_length[ref]),
                "interval": int(interval),
                "values": []
            }

            # array to store the values of coverage for each interval
            array_of_cov = []
            # the counter that is used to output the values per interval
            reset_counter = 0
            # loop to generate dict_cov
            logger.info("Generating plot data for plasmid: {}".format(ref))
            for i in range(int(plasmid_length[ref])):
                # checks if key for a given position is in dict and if so
                # adds it to array of cov, otherwise it will add a 0
                try:
                    array_of_cov.append(int(depth_dic_coverage[ref][str(i)]))
                except KeyError:
                    array_of_cov.append(0)

                # if the counter equals the interval then output to dict_cov
                if reset_counter == interval:
                    dict_cov[ref]["values"].append(
                        int(sum(array_of_cov)/len(array_of_cov))
                    )
                    # reset counter
                    reset_counter = 0
                else:
                    # if counter is less than interval then sums 1
                    reset_counter += 1

    logger.info("Successfully generated dicts necessary for output json file "
                "and .report.json depth file.")
    logger.debug("Size of percentage_bases_covered: {} kb".format(
        asizeof(percentage_bases_covered)/1024))
    logger.debug("Size of dict_cov: {} kb".format(asizeof(dict_cov)/1024))
    return percentage_bases_covered, dict_cov


@MainWrapper
def main(depth_file, json_dict, cutoff, sample_id):
    """
    Function that handles the inputs required to parse depth files from bowtie
    and dumps a dict to a json file that can be imported into pATLAS.

    Parameters
    ----------
    depth_file: str
         the path to depth file for each sample
    json_dict: str
        the file that contains the dictionary with keys and values for
        accessions
        and their respective lengths
    cutoff: str
        the cutoff used to trim the unwanted matches for the minimum coverage
        results from mapping. This value may range between 0 and 1.
    sample_id: str
        the id of the sample being parsed

    """

    # check for the appropriate value for the cutoff value for coverage results
    logger.debug("Cutoff value: {}. Type: {}".format(cutoff, type(cutoff)))
    try:
        cutoff_val = float(cutoff)
        if cutoff_val < 0.4:
            logger.warning("This cutoff value will generate a high volume of "
                           "plot data. Therefore '.report.json' can be too big")
    except ValueError:
        logger.error("Cutoff value should be a string such as: '0.6'. "
                     "The outputted value: {}. Make sure to provide an "
                     "appropriate value for --cov_cutoff".format(cutoff))
        sys.exit(1)

    # loads dict from file, this file is provided in docker image

    plasmid_length = json.load(open(json_dict))
    if plasmid_length:
        logger.info("Loaded dictionary of plasmid lengths")
    else:
        logger.error("Something went wrong and plasmid lengths dictionary"
                     "could not be loaded. Check if process received this"
                     "param successfully.")
        sys.exit(1)

    # read depth file
    depth_file_in = open(depth_file)

    # first reads the depth file and generates dictionaries to handle the input
    # to a simpler format
    logger.info("Reading depth file and creating dictionary to dump.")
    depth_dic_coverage = depth_file_reader(depth_file_in)
    percentage_bases_covered, dict_cov = generate_jsons(depth_dic_coverage,
                                                        plasmid_length,
                                                        cutoff_val)

    if percentage_bases_covered and dict_cov:
        logger.info("percentage_bases_covered length: {}".format(
            str(len(percentage_bases_covered))))
        logger.info("dict_cov length: {}".format(str(len(dict_cov))))
    else:
        logger.error("Both dicts that dump to JSON file or .report.json are "
                     "empty.")

    # then dump do file
    logger.info("Dumping to {}".format("{}_mapping.json".format(depth_file)))
    with open("{}_mapping.json".format(depth_file), "w") as output_json:
        output_json.write(json.dumps(percentage_bases_covered))

    json_dic = {
        "tableRow": [{
            "sample": sample_id,
            "data": [{
                "header": "Mapping",
                "table": "plasmids",
                "patlas_mapping": percentage_bases_covered,
                "value": len(percentage_bases_covered)
            }]
        }],
        "sample": sample_id,
        "patlas_mapping": percentage_bases_covered,
        "plotData": [{
            "sample": sample_id,
            "data": {
                "patlasMappingSliding": dict_cov
            },
        }]
    }

    logger.debug("Size of dict_cov: {} kb".format(asizeof(json_dic)/1024))
    logger.info("Writing to .report.json")
    with open(".report.json", "w") as json_report:
        json_report.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == "__main__":
    main(DEPTH_TXT, JSON_LENGTH, CUTOFF, SAMPLE_ID)
