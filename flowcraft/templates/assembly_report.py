#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to provide a summary report for a given assembly
in Fasta format.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``assembly`` : Path to assembly file in Fasta format.
    - e.g.: ``'assembly.fasta'``

Generated output
----------------

- ``${sample_id}_assembly_report.csv`` : CSV with summary information of the \
    assembly.
    - e.g.: ``'SampleA_assembly_report.csv'``

Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "16012018"
__template__ = "assembly_report-nf"

import os
import re
import json
import traceback
import subprocess

from collections import OrderedDict
from subprocess import PIPE

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


def __get_version_pilon():

    pilon_path = "/NGStools/pilon-1.22.jar"

    try:

        cli = ["java", "-jar", pilon_path, "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        version = stdout.split()[2].decode("utf8")

    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return {
        "program": "Pilon",
        "version": version,
    }


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLY_FILE = '$assembly'
    COVERAGE_BP_FILE = '$coverage_bp'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLY_FILE: {}".format(ASSEMBLY_FILE))
    logger.debug("COVERAGE_BP_FILE: {}".format(COVERAGE_BP_FILE))


class Assembly:
    """Class that parses and filters an assembly file in Fasta format.

    This class parses an assembly file, collects a number
    of summary statistics and metadata from the contigs and reports.

    Parameters
    ----------
    assembly_file : str
        Path to assembly file.
    sample_id : str
        Name of the sample for the current assembly.
    """

    def __init__(self, assembly_file, sample_id):

        self.summary_info = OrderedDict([
            ("ncontigs", 0),
            ("avg_contig_size", []),
            ("n50", 0),
            ("total_len", 0),
            ("avg_gc", []),
            ("missing_data", 0)
        ])
        """
        OrderedDict: Initialize summary information dictionary. Contains keys:

            - ``ncontigs``: Number of contigs
            - ``avg_contig_size``: Average size of contigs
            - ``n50``: N50 metric
            - ``total_len``: Total assembly length
            - ``avg_gc``: Average GC proportion
            - ``missing_data``: Count of missing data characters
        """

        self.contigs = OrderedDict()
        """
        OrderedDict: Object that maps the contig headers to the corresponding
        sequence
        """

        self.contig_coverage = OrderedDict()
        """
        OrderedDict: Object that maps the contig headers to the corresponding
        list of per-base coverage
        """

        self.sample = sample_id
        """
        str: Sample id
        """

        self.contig_boundaries = {}
        """
        dict: Maps the boundaries of each contig in the genome
        """

        self._parse_assembly(assembly_file)

    def _parse_assembly(self, assembly_file):
        """Parse an assembly file in fasta format.

        This is a Fasta parsing method that populates the
        :py:attr:`Assembly.contigs` attribute with data for each contig in the
         assembly.

        Parameters
        ----------
        assembly_file : str
            Path to the assembly fasta file.

        """

        with open(assembly_file) as fh:

            header = None
            logger.debug("Starting iteration of assembly file: {}".format(
                assembly_file))

            for line in fh:

                # Skip empty lines
                if not line.strip():
                    continue

                if line.startswith(">"):
                    # Add contig header to contig dictionary
                    header = line[1:].strip()
                    self.contigs[header] = []

                else:
                    # Add sequence string for the current contig
                    self.contigs[header].append(line.strip())

            # After populating the contigs dictionary, convert the values
            # list into a string sequence
            self.contigs = OrderedDict(
                (header, "".join(seq)) for header, seq in self.contigs.items())

    @staticmethod
    def _get_contig_id(contig_str):
        """Tries to retrieve contig id. Returns the original string if it
        is unable to retrieve the id.

        Parameters
        ----------
        contig_str : str
            Full contig string (fasta header)

        Returns
        -------
        str
            Contig id
        """

        contig_id = contig_str

        try:
            contig_id = re.search(".*NODE_([0-9]*)_.*", contig_str).group(1)
        except AttributeError:
            pass

        try:
            contig_id = re.search(".*Contig_([0-9]*)_.*", contig_str).group(1)
        except AttributeError:
            pass

        return contig_id

    def get_summary_stats(self, output_csv=None):
        """Generates a CSV report with summary statistics about the assembly

        The calculated statistics are:

            - Number of contigs
            - Average contig size
            - N50
            - Total assembly length
            - Average GC content
            - Amount of missing data

        Parameters
        ----------
        output_csv: str
            Name of the output CSV file.
        """

        contig_size_list = []

        self.summary_info["ncontigs"] = len(self.contigs)

        for contig_id, sequence in self.contigs.items():

            logger.debug("Processing contig: {}".format(contig_id))

            # Get contig sequence size
            contig_len = len(sequence)

            # Add size for average contig size
            contig_size_list.append(contig_len)

            # Add to total assembly length
            self.summary_info["total_len"] += contig_len

            # Add to average gc
            self.summary_info["avg_gc"].append(
                sum(map(sequence.count, ["G", "C"])) / contig_len
            )

            # Add to missing data
            self.summary_info["missing_data"] += sequence.count("N")

        # Get average contig size
        logger.debug("Getting average contig size")
        self.summary_info["avg_contig_size"] = \
            sum(contig_size_list) / len(contig_size_list)

        # Get average gc content
        logger.debug("Getting average GC content")
        self.summary_info["avg_gc"] = \
            sum(self.summary_info["avg_gc"]) / len(self.summary_info["avg_gc"])

        # Get N50
        logger.debug("Getting N50")
        cum_size = 0
        for l in sorted(contig_size_list, reverse=True):
            cum_size += l
            if cum_size >= self.summary_info["total_len"] / 2:
                self.summary_info["n50"] = l
                break

        if output_csv:
            logger.debug("Writing report to csv")
            # Write summary info to CSV
            with open(output_csv, "w") as fh:
                summary_line = "{}, {}\\n".format(
                    self.sample, ",".join(
                        [str(x) for x in self.summary_info.values()]))
                fh.write(summary_line)

    def _get_window_labels(self, window):
        """Returns the mapping between sliding window points and their contigs,
        and the x-axis position of contig

        Parameters
        ----------
        window : int
            Size of the window.

        Returns
        -------
        xbars : list
            The x-axis position of the ending for each contig.
        labels : list
            The x-axis labels for each data point in the sliding window

        """

        # Get summary stats, if they have not yet been triggered
        if not self.summary_info:
            self.get_summary_stats()

        # Get contig boundary positon
        c = 0
        xbars = []
        for contig, seq in self.contigs.items():
            contig_id = self._get_contig_id(contig)
            self.contig_boundaries[contig_id] = [c, c + len(seq)]
            c += len(seq)
            xbars.append((contig_id, c, contig))

        return xbars

    @staticmethod
    def _gc_prop(s, length):
        """Get proportion of GC from a string

        Parameters
        ----------
        s : str
            Arbitrary string

        Returns
        -------
        x : float
            GC proportion.
        """

        gc = sum(map(s.count, ["c", "g"]))

        return gc / length

    def get_gc_sliding(self, window=2000):
        """Calculates a sliding window of the GC content for the assembly


        Returns
        -------
        gc_res : list
            List of GC proportion floats for each data point in the sliding
            window
        """

        gc_res = []

        # Get complete sequence to calculate sliding window values
        complete_seq = "".join(self.contigs.values()).lower()

        for i in range(0, len(complete_seq), window):

            seq_window = complete_seq[i:i + window]

            # Get GC proportion
            gc_res.append(round(self._gc_prop(seq_window, len(seq_window)), 2))

        return gc_res

    def _get_coverage_from_file(self, coverage_file):
        """

        Parameters
        ----------
        coverage_file

        Returns
        -------

        """

        with open(coverage_file) as fh:

            for line in fh:

                fields = line.strip().split()

                # Get header
                header = fields[0]
                coverage = int(fields[2])

                if header not in self.contig_coverage:
                    self.contig_coverage[header] = [coverage]
                else:
                    self.contig_coverage[header].append(coverage)

    def get_coverage_sliding(self, coverage_file, window=2000):
        """

        Parameters
        ----------
        coverage_file : str
            Path to file containing the coverage info at the per-base level
            (as generated by samtools depth)
        window : int
            Size of sliding window

        Returns
        -------

        """

        if not self.contig_coverage:
            self._get_coverage_from_file(coverage_file)

        # Stores the coverage results
        cov_res = []

        # Make flat list of coverage values across genome
        complete_cov = [x for y in self.contig_coverage.values() for x in y]

        for i in range(0, len(complete_cov), window):
            # Get coverage values for current window
            cov_window = complete_cov[i:i + window]
            # Get mean coverage
            cov_res.append(int(sum(cov_window) / len(cov_window)))

        return cov_res


@MainWrapper
def main(sample_id, assembly_file, coverage_bp_file=None):
    """Main executor of the assembly_report template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    assembly_file : str
        Path to assembly file in Fasta format.

    """

    logger.info("Starting assembly report")
    assembly_obj = Assembly(assembly_file, sample_id)

    logger.info("Retrieving summary statistics for assembly")
    assembly_obj.get_summary_stats("{}_assembly_report.csv".format(sample_id))

    size_dist = [len(x) for x in assembly_obj.contigs.values()]
    json_dic = {
        "tableRow": [{
            "sample": sample_id,
            "data": [
                {"header": "Contigs",
                 "value": assembly_obj.summary_info["ncontigs"],
                 "table": "assembly",
                 "columnBar": True},
                {"header": "Assembled BP",
                 "value": assembly_obj.summary_info["total_len"],
                 "table": "assembly",
                 "columnBar": True},
            ]
        }],
        "plotData": [{
            "sample": sample_id,
            "data": {
                "size_dist": size_dist
            }
        }]
    }

    if coverage_bp_file:
        try:
            window = 2000
            gc_sliding_data = assembly_obj.get_gc_sliding(window=window)
            cov_sliding_data = \
                assembly_obj.get_coverage_sliding(coverage_bp_file,
                                                  window=window)

            # Get total basepairs based on the individual coverage of each
            # contig bpx
            total_bp = sum(
                [sum(x) for x in assembly_obj.contig_coverage.values()]
            )

            # Add data to json report
            json_dic["plotData"][0]["data"]["genomeSliding"] = {
                "gcData": gc_sliding_data,
                "covData": cov_sliding_data,
                "window": window,
                "xbars": assembly_obj._get_window_labels(window),
                "assemblyFile": os.path.basename(assembly_file)
            }
            json_dic["plotData"][0]["data"]["sparkline"] = total_bp

        except:
            logger.error("Unexpected error creating sliding window data:\\n"
                         "{}".format(traceback.format_exc()))

    # Write json report
    with open(".report.json", "w") as json_report:

        json_report.write(json.dumps(json_dic, separators=(",", ":")))

    with open(".status", "w") as status_fh:
        status_fh.write("pass")


if __name__ == '__main__':

    main(SAMPLE_ID, ASSEMBLY_FILE, COVERAGE_BP_FILE)

