#!/usr/bin/env python3

"""
Purpose
-------

This module is intended parse the results of the Trimmomatic log for a set
of one or more samples.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``log_files``: Trimmomatic log files.
    - e.g.: ``'Sample1_trimlog.txt Sample2_trimlog.txt'``


Generated output
----------------
- ``trimmomatic_report.csv`` : Summary report of the trimmomatic logs for\
    all samples

Code documentation
------------------

"""

__version__ = "1.0.0"
__build__ = "16012018"
__template__ = "trimmomatic_report-nf"

import os
import json

from collections import OrderedDict

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


if __file__.endswith(".command.sh"):
    LOG_FILES = '$log_files'.split()


def parse_log(log_file):
    """Retrieves some statistics from a single Trimmomatic log file.

    This function parses Trimmomatic's log file and stores some trimming
    statistics in an :py:class:`OrderedDict` object. This object contains
    the following keys:

        - ``clean_len``: Total length after trimming.
        - ``total_trim``: Total trimmed base pairs.
        - ``total_trim_perc``: Total trimmed base pairs in percentage.
        - ``5trim``: Total base pairs trimmed at 5' end.
        - ``3trim``: Total base pairs trimmed at 3' end.

    Parameters
    ----------
    log_file : str
        Path to trimmomatic log file.

    Returns
    -------
    x : :py:class:`OrderedDict`
        Object storing the trimming statistics.

    """

    template = OrderedDict([
        # Total length after trimming
        ("clean_len", 0),
        # Total trimmed base pairs
        ("total_trim", 0),
        # Total trimmed base pairs in percentage
        ("total_trim_perc", 0),
        # Total trimmed at 5' end
        ("5trim", 0),
        # Total trimmed at 3' end
        ("3trim", 0),
        # Bad reads (completely trimmed)
        ("bad_reads", 0)
    ])

    with open(log_file) as fh:

        for line in fh:
            # This will split the log fields into:
            # 0. read length after trimming
            # 1. amount trimmed from the start
            # 2. last surviving base
            # 3. amount trimmed from the end
            fields = [int(x) for x in line.strip().split()[-4:]]

            if not fields[0]:
                template["bad_reads"] += 1

            template["5trim"] += fields[1]
            template["3trim"] += fields[3]
            template["total_trim"] += fields[1] + fields[3]
            template["clean_len"] += fields[0]

        total_len = template["clean_len"] + template["total_trim"]

        if total_len:
            template["total_trim_perc"] = round(
                (template["total_trim"] / total_len) * 100, 2)
        else:
            template["total_trim_perc"] = 0

    return template


def write_report(storage_dic, output_file, sample_id):
    """ Writes a report from multiple samples.

    Parameters
    ----------
    storage_dic : dict or :py:class:`OrderedDict`
        Storage containing the trimming statistics. See :py:func:`parse_log`
        for its generation.
    output_file : str
        Path where the output file will be generated.
    sample_id : str
        Id or name of the current sample.
    """

    with open(output_file, "w") as fh, open(".report.json", "w") as json_rep:

        # Write header
        fh.write("Sample,Total length,Total trimmed,%,5end Trim,3end Trim,"
                 "bad_reads\\n")

        # Write contents
        for sample, vals in storage_dic.items():
            fh.write("{},{}\\n".format(
                sample, ",".join([str(x) for x in vals.values()])))

            json_dic = {
                "tableRow": [{
                    "sample": sample_id,
                    "data": [
                        {"header": "trimmed",
                         "value": vals["total_trim_perc"],
                         "table": "qc",
                         "columnBar": True},
                    ]
                }],
                "plotData": [{
                    "sample": sample_id,
                    "data": {
                        "sparkline": vals["clean_len"]
                    }
                }],
                "badReads": vals["bad_reads"]
            }
            json_rep.write(json.dumps(json_dic, separators=(",", ":")))


@MainWrapper
def main(log_files):
    """ Main executor of the trimmomatic_report template.

    Parameters
    ----------
    log_files : list
        List of paths to the trimmomatic log files.
    """

    log_storage = OrderedDict()

    for log in log_files:

        log_id = log.rstrip("_trimlog.txt")

        # Populate storage of current sample
        log_storage[log_id] = parse_log(log)

        # Remove temporary trim log file
        os.remove(log)

    write_report(log_storage, "trimmomatic_report.csv", log_id)


if __name__ == '__main__':

    main(LOG_FILES)
