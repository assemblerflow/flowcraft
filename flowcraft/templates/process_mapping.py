#!/usr/bin/env python3

import re
import os
import json

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper



"""
Purpose
-------

This module is intended to process the output of mapping proces from a single
sample from the program Bowtie for the report component.
The main input is an log file produced by the mapper.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id``: Sample Identification string.
    - e.g.: ``'SampleA'``
- ``bowtie_log``: Log file from the mapper.
    - e.g.: ``'bowtie.log'``

Generated output
----------------
- ``.report.jason``: Data structure for the report

Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "10.09.2018"
__template__ = "remove_host-nf"

logger = get_logger(__file__)


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    BOWTIE_LOG = '$bowtie_log'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("BOWTIE_LOG: {}".format(BOWTIE_LOG))



class Bowtie:
    """
    Class to parse and store the info in the bowtie log file.

    """

    def __init__(self, sample_id, bowtie_log):

        self.sample = sample_id
        """
        str: The name of the sample for the assembly.
        """

        self.n_reads = 0

        self.align_0x = 0

        self.align_1x = 0

        self.align_mt1x = 0

        self.overall_rate = 0.0

        # Parse assembly and populate self.n_reads, self.align_0x, self.align_1x, self.align_mt1x and self.overall_rate
        self.parse_log(bowtie_log)


    def set_n_reads(self, n_reads):
        self.n_reads = int(n_reads)


    def set_align_0x(self,align_0x):
        self.align_0x = align_0x


    def set_align_1x(self,align_1x):
        self.align_1x = align_1x


    def set_align_mt1x(self,align_mt1x):
        self.align_mt1x = align_mt1x


    def set_overall_rate(self,overall_rate):
        self.overall_rate = overall_rate


    def parse_log(self, bowtie_log):
        """Parse a bowtie log file.

        This is a bowtie log parsing method that populates the
        :py:attr:`self.n_reads, self.align_0x, self.align_1x, self.align_mt1x and self.overall_rate` attributes with
        data from the log file.

        Disclamer: THIS METHOD IS HORRIBLE BECAUSE THE BOWTIE LOG IS HORRIBLE.

        The insertion of data on the attribytes is done by the
        :py:meth:`set_attribute method.

        Parameters
        ----------
        bowtie_log : str
            Path to the boetie log file.

       """

        print("is here!")

        # Regexes - thanks to https://github.com/ewels/MultiQC/blob/master/multiqc/modules/bowtie2/bowtie2.py
        regexes = {
            'unpaired': {
                'unpaired_aligned_none': r"(\\d+) \\([\\d\\.]+%\\) aligned 0 times",
                'unpaired_aligned_one': r"(\\d+) \\([\\d\\.]+%\\) aligned exactly 1 time",
                'unpaired_aligned_multi': r"(\\d+) \\([\\d\\.]+%\\) aligned >1 times"
            },
            'paired': {
                'paired_aligned_none': r"(\\d+) \\([\\d\\.]+%\\) aligned concordantly 0 times",
                'paired_aligned_one': r"(\\d+) \\([\\d\\.]+%\\) aligned concordantly exactly 1 time",
                'paired_aligned_multi': r"(\\d+) \\([\\d\\.]+%\\) aligned concordantly >1 times",
                'paired_aligned_discord_one': r"(\\d+) \\([\\d\\.]+%\\) aligned discordantly 1 time",
                'paired_aligned_discord_multi': r"(\\d+) \\([\\d\\.]+%\\) aligned discordantly >1 times",
                'paired_aligned_mate_one': r"(\\d+) \\([\\d\\.]+%\\) aligned exactly 1 time",
                'paired_aligned_mate_multi': r"(\\d+) \\([\\d\\.]+%\\) aligned >1 times",
                'paired_aligned_mate_none': r"(\\d+) \\([\\d\\.]+%\\) aligned 0 times"
            }
        }

        #Missing parser for unpaired (not implemented in flowcraft yet)

        with open(bowtie_log, "r") as f:
            #Go through log file line by line
            for l in f:

                print(l)

                #total reads
                total = re.search(r"(\\d+) reads; of these:", l)
                print(total)
                if total:
                    print(total)
                    self.set_n_reads(total.group(1))


                # Paired end reads aka the pain
                paired = re.search(r"(\\d+) \\([\\d\\.]+%\\) were paired; of these:", l)
                if paired:
                    paired_total = int(paired.group(1))

                    paired_numbers = {}

                    # Do nested loop whilst we have this level of indentation
                    l = f.readline()
                    while l.startswith('    '):
                        for k, r in regexes['paired'].items():
                            match = re.search(r, l)
                            if match:
                                paired_numbers[k] = int(match.group(1))
                        l = f.readline()


                    align_zero_times = paired_numbers['paired_aligned_none'] + paired_numbers['paired_aligned_mate_none']
                    if align_zero_times:
                        self.set_align_0x(align_zero_times)

                    align_one_time = paired_numbers['paired_aligned_one'] + paired_numbers['paired_aligned_mate_one']
                    if align_one_time:
                        self.set_align_1x(align_one_time)

                    align_more_than_one_time = paired_numbers['paired_aligned_multi'] + paired_numbers['paired_aligned_mate_multi']
                    if align_more_than_one_time:
                        self.set_align_mt1x(align_more_than_one_time)


                # Overall alignment rate
                overall = re.search(r"([\\d\\.]+)% overall alignment rate", l)
                if overall:
                    self.overall_rate = float(overall.group(1))


@MainWrapper
def main(sample_id, bowite_log):
    """Main executor of the process_mapping template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    boetie_log: str
        Path to the log file generated by bowtie.

    """

    logger.info("Starting mapping file processing")
    warnings = []
    fails = ""

    bowtie_info = Bowtie(sample_id, bowite_log)

    print(bowtie_info.overall_rate)


    with open(".report.json", "w") as json_report:
        json_dic = {
            "tableRow": [{
                "sample": sample_id,
                "data": [
                    {"header": "Reads",
                     "value": int(bowtie_info.n_reads),
                     "table": "mapping",
                     "columnBar": False},
                    {"header": "Unmapped",
                     "value": int(bowtie_info.align_0x),
                     "table": "mapping",
                     "columnBar": False},
                    {"header": "Mapped 1x",
                     "value": int(bowtie_info.align_1x),
                     "table": "mapping",
                     "columnBar": False},
                    {"header": "Mapped >1x",
                     "value": int(bowtie_info.align_mt1x),
                     "table": "mapping",
                     "columnBar": False},
                    {"header": "Overall alignment rate (%)",
                     "value": float(bowtie_info.overall_rate),
                     "table": "mapping",
                     "columnBar": False}
                ]
            }],
        }

        if warnings:
            json_dic["warnings"] = [{
                "sample": sample_id,
                "table": "mapping",
                "value": warnings
            }]

        if fails:
            json_dic["fail"] = [{
                "sample": sample_id,
                "table": "mapping",
                "value": [fails]
            }]

        json_report.write(json.dumps(json_dic, separators=(",", ":")))

    with open(".status", "w") as status_fh:
        status_fh.write("pass")


if __name__ == '__main__':

    main(SAMPLE_ID, BOWTIE_LOG)