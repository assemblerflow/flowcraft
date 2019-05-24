#!/usr/bin/env python3

"""
Purpose
-------
This module is intended to process the output of concoct
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
    CLUSTER = '$cluster'
    CONTIGS = '$contigs'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("CLUSTER: {}".format(CLUSTER))
    logger.debug("CONTIGS: {}".format(CONTIGS))


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


def parse_cluster_csv(file):
    """
    Simple csv parser for clustering file of concoct
    :param file: clustering csv file
    :return: dictionary containing the cluster id and the contigs in the cluster
    """

    clusters = {}

    reader = csv.reader(open(file), delimiter=',')
    next(reader)  # skip header
    for row in reader:
        if row[1] in clusters:
            clusters[row[1]].append(row[0])
        else:
            clusters[row[1]] = [row[0]]

    return clusters


def get_GC(sequence):

    return round(sum(1 for nucl in sequence if nucl in ['G', 'C'])/len(sequence)*100, 2)


def merge_data(contigs, clusters):
    """
    Obtain genome size, cg content and number of contigs for concoct bins

    :param contigs: dict with the sequences for the binned contigs
    :param clusters: dict with the cluster and respective sequence headers
    :return: dict with the statistics for each bin (cluster)
    """

    binning = {}

    for cluster_id in clusters.keys():
        complete_sequence = ''
        n_sequences = 0
        for sequence in clusters[cluster_id]:
            complete_sequence += contigs[sequence]
            n_sequences += 1

        binning[int(cluster_id)] = {"Bin name": cluster_id,
                                    "Contig number": n_sequences,
                                    "Genome size": len(complete_sequence),
                                    "GC content": get_GC(complete_sequence)}

    return binning



@MainWrapper
def main(sample_id, cluster_file, contig_file):

    seqs = parse_assembly(contig_file)

    clusters = parse_cluster_csv(cluster_file)

    bin_stats = merge_data(seqs, clusters)

    report_list = [["Bin name", "Contig number", "Genome size", "GC content %"]]

    for key, value in sorted(bin_stats.items(), key=lambda x: x[0]):
        print("{} : {}".format(key, value))
        report_list.append([value["Bin name"],
                            str(value["Contig number"]),
                            str(value["Genome size"]),
                            str(value["GC content"])])

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
    main(SAMPLE_ID, CLUSTER, CONTIGS)
