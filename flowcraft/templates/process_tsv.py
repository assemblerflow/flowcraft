#!/usr/bin/env python3

"""
Purpose
-------
This module is intended to process the output in tsv
 to generate a report in json format.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``sample_id`` : Sample Identification string.
- ``tsv``: tsv output.

"""

import json
import csv
import os

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

__version__ = "1.0.1"
__build__ = "05.10.2018"
__template__ = "maxbin2-nf"

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    FILE = '$tsv'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("FILE: {}".format(FILE))

@MainWrapper
def main(sample_id, tsv_file):

    # this tsvData could be a single object since it only has one element
    # this data type expects full tables in tsv format
    report_json = {
        "tsvData": [{
            "sample": sample_id,
            "data": {}
        }]
    }

    # web-app excepts a list with all the values in the table.
    #  To expand this to other processes other than MaxBin2, this line needs to be reworked
    report_json["tsvData"][0]["data"]["MaxBin2"] = list(csv.reader(open(tsv_file), delimiter='\t'))

    with open(".report.json", "w") as k:
        k.write(json.dumps(report_json))


if __name__ == "__main__":
    main(SAMPLE_ID, FILE)
