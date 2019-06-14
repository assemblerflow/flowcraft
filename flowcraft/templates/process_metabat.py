#!/usr/bin/env python3

"""
Purpose
-------
This module is intended to process the output of metaBAT
 to generate a report in json format.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``sample_id`` : Sample Identification string.
- ``cluster``: concoct cluster output.

"""

import json
import csv
import os
from itertools import groupby

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper


__version__ = "1.0.0"
__build__ = "22.05.2019"
__template__ = "concoct-nf"

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    BINS = '$bins'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("BINS: {}".format(BINS))


def parse_assembly(file):
    """
    Simple fasta parser.
    :param file: assembly file in fasta format
    :return: dictionary containing the contigs in the assembly
    """

    all_seqs = {}

    with open(file, "r") as handle:
        entry = (x[1] for x in groupby(handle, lambda line: line[0] == ">"))
        for header in entry:
            contig_header = header.__next__()[1:].strip()
            contig_seq = "".join(s.strip() for s in entry.__next__())
            all_seqs[contig_header] = contig_seq

    return all_seqs


def get_cg(sequence):

    return round(sum(1 for nucl in sequence if nucl in ['G', 'C'])/len(sequence)*100, 2)


def get_bin_stats(bin_file):
    n_contigs = 0
    all_seq = ""

    with open(bin_file, "r") as handle:
        entry = (x[1] for x in groupby(handle, lambda line: line[0] == ">"))
        for header in entry:
            n_contigs += 1
            all_seq += "".join(s.strip() for s in entry.__next__())

    return str(n_contigs), str(len(all_seq)), str(get_cg(all_seq))

@MainWrapper
def main(sample_id, bins):

    report_list = [["Bin name", "Contig number", "Genome size", "GC content %"]]

    if len(bins) == 1 and "false_bin.fa" not in bins:
        ncontigs, gsize, gc = get_bin_stats(bins)
        report_list.append([bins.split(".")[1], ncontigs, gsize, gc])
    else:
        for file in bins:
            ncontigs, gsize, gc = get_bin_stats(file)
            report_list.append([file.split(".")[1], ncontigs, gsize, gc])

    # this tsvData is a single object since it only has one element
    # this data type expects full tables in tsv format
    report_json = {
        "tsvData": [{
            "sample": sample_id,
            "data": {}
        }]
    }

    # web-app excepts a list with all the values in the table.
    #  To expand this to other processes other than MaxBin2, this line needs to be reworked
    report_json["tsvData"][0]["data"]["MaxBin2"] = report_list

    with open(".report.json", "w") as k:
        k.write(json.dumps(report_json))


if __name__ == "__main__":
    main(SAMPLE_ID, BINS)
