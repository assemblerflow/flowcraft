#!/usr/bin/env python3

"""
Purpose
-------

This module is intended parse the results of FastQC for paired end FastQ \
samples. It parses two reports:

    - Categorical report
    - Nucleotide level report.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample identification string
    - e.g.: ``'SampleA'``

- ``result_p1`` : Path to both FastQC result files for pair 1
    - e.g.: ``'SampleA_1_data SampleA_1_summary'``

- ``result_p2`` : Path to both FastQC result files for pair 2
    - e.g.: ``'SampleA_2_data SampleA_2_summary'``

- ``opts`` : *Specify additional arguments for executing fastqc_report. \
    The arguments should be a string of command line arguments,\
    The accepted arguments are:*
    - ``'--ignore-tests'`` : Ignores test results from FastQC categorical\
    summary. This is used in the first run of FastQC.

Generated output
----------------

The generated output are output files that contain an object, usually a string.

- ``fastqc_health`` : Stores the health check for the current sample. If it
    passes all checks, it contains only the string 'pass'. Otherwise, contains
    the summary categories and their respective results
    - e.g.: ``'pass'``
- ``optimal_trim`` : Stores a tuple with the optimal trimming positions for 5'
    and 3' ends of the reads.
    - e.g.: ``'15 151'``

Code documentation
------------------

"""

__version__ = "1.0.2"
__build__ = "12052018"
__template__ = "fastqc_report-nf"

import os
import json

from collections import OrderedDict

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    RESULT_P1 = '$result_p1'.split()
    RESULT_P2 = '$result_p2'.split()
    SAMPLE_ID = '$sample_id'
    OPTS = '$opts'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("RESULT_P1: {}".format(RESULT_P1))
    logger.debug("RESULT_P2: {}".format(RESULT_P2))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("OPTS: {}".format(OPTS))


def _get_quality_stats(d, start_str, field_start=1, field_end=2):
    """

    Parameters
    ----------
    d

    Returns
    -------

    """

    min_parsed = False
    parse = False
    report = []
    start_str = start_str
    end_str = ">>END_MODULE"

    with open(d) as fh:

        for line in fh:

            if line.startswith(start_str):
                next(fh)
                parse = True
                status = line.strip().split()[-1]

            # Exit parser when end string is found
            elif parse and line.startswith(end_str):
                return report, status

            elif parse:

                fields = line.strip().split()

                # This is triggered when the first value of a line series is
                # not 1. If the starting point of the series is a number
                # different from 1, fill the report with 0 until that point
                if not min_parsed:
                    if fields[0] != "1":
                        try:
                            blank_points = int(fields[0]) - 1
                            report.extend([0] * blank_points)
                        except ValueError:
                            pass
                    min_parsed = True

                report.append(";".join([
                    str(round(float(x), 2)) for x in
                    fields[field_start: field_end]
                ]))


def write_json_report(sample_id, data1, data2):
    """Writes the report

    Parameters
    ----------
    data1
    data2

    Returns
    -------

    """

    parser_map = {
        "base_sequence_quality": ">>Per base sequence quality",
        "sequence_quality": ">>Per sequence quality scores",
        "base_gc_content": ">>Per sequence GC content",
        "base_n_content": ">>Per base N content",
        "sequence_length_dist": ">>Sequence Length Distribution",
        "per_base_sequence_content": ">>Per base sequence content"
    }

    json_dic = {
        "plotData": [{
            "sample": sample_id,
            "data": {
                "base_sequence_quality": {"status": None, "data": []},
                "sequence_quality": {"status": None, "data": []},
                "base_gc_content": {"status": None, "data": []},
                "base_n_content": {"status": None, "data": []},
                "sequence_length_dist": {"status": None, "data": []},
                "per_base_sequence_content": {"status": None, "data": []}
            }
        }]
    }

    for cat, start_str in parser_map.items():

        if cat == "per_base_sequence_content":
            fs = 1
            fe = 5
        else:
            fs = 1
            fe = 2

        report1, status1 = _get_quality_stats(data1, start_str,
                                              field_start=fs, field_end=fe)
        report2, status2 = _get_quality_stats(data2, start_str,
                                              field_start=fs, field_end=fe)

        status = None
        for i in ["fail", "warn", "pass"]:
            if i in [status1, status2]:
                status = i

        json_dic["plotData"][0]["data"][cat]["data"] = [report1, report2]
        json_dic["plotData"][0]["data"][cat]["status"] = status

    return json_dic


def get_trim_index(biased_list):
    """Returns the trim index from a ``bool`` list

    Provided with a list of ``bool`` elements (``[False, False, True, True]``),
    this function will assess the index of the list that minimizes the number
    of True elements (biased positions) at the extremities. To do so,
    it will iterate over the boolean list and find an index position where
    there are two consecutive ``False`` elements after a ``True`` element. This
    will be considered as an optimal trim position. For example, in the
    following list::

        [True, True, False, True, True, False, False, False, False, ...]

    The optimal trim index will be the 4th position, since it is the first
    occurrence of a ``True`` element with two False elements after it.

    If the provided ``bool`` list has no ``True`` elements, then the 0 index is
    returned.

    Parameters
    ----------
    biased_list: list
        List of ``bool`` elements, where ``True`` means a biased site.

    Returns
    -------
        x : index position of the biased list for the optimal trim.

    """

    # Return index 0 if there are no biased positions
    if set(biased_list) == {False}:
        return 0

    if set(biased_list[:5]) == {False}:
        return 0

    # Iterate over the biased_list array. Keep the iteration going until
    # we find a biased position with the two following positions unbiased
    # (e.g.: True, False, False).
    # When this condition is verified, return the last biased position
    # index for subsequent trimming.
    for i, val in enumerate(biased_list):
        if val and set(biased_list[i+1:i+3]) == {False}:
            return i + 1

    # If the previous iteration could not find and index to trim, it means
    # that the whole list is basically biased. Return the length of the
    # biased_list
    return len(biased_list)


def trim_range(data_file):
    """Assess the optimal trim range for a given FastQC data file.

    This function will parse a single FastQC data file, namely the
    *'Per base sequence content'* category. It will retrieve the A/T and G/C
    content for each nucleotide position in the reads, and check whether the
    G/C and A/T proportions are between 80% and 120%. If they are, that
    nucleotide position is marked as biased for future removal.

    Parameters
    ----------
    data_file: str
        Path to FastQC data file.

    Returns
    -------
    trim_nt: list
        List containing the range with the best trimming positions for the
        corresponding FastQ file. The first element is the 5' end trim index
        and the second element is the 3' end trim index.
    """

    logger.debug("Starting trim range assessment")

    # Target string for nucleotide bias assessment
    target_nuc_bias = ">>Per base sequence content"
    logger.debug("Target string to start nucleotide bias assessment set to "
                 "{}".format(target_nuc_bias))
    # This flag will become True when gathering base proportion data
    # from file.
    gather = False

    # This variable will store a boolean array on the biased/unbiased
    # positions. Biased position will be True, while unbiased positions
    # will be False
    biased = []

    with open(data_file) as fh:

        for line in fh:
            # Start assessment of nucleotide bias
            if line.startswith(target_nuc_bias):
                # Skip comment line
                logger.debug("Found target string at line: {}".format(line))
                next(fh)
                gather = True
            # Stop assessment when reaching end of target module
            elif line.startswith(">>END_MODULE") and gather:
                logger.debug("Stopping parsing at line: {}".format(line))
                break
            elif gather:
                # Get proportions of each nucleotide
                g, a, t, c = [float(x) for x in line.strip().split()[1:]]
                # Get 'GC' and 'AT content
                gc = (g + 0.1) / (c + 0.1)
                at = (a + 0.1) / (t + 0.1)
                # Assess bias
                if 0.8 <= gc <= 1.2 and 0.8 <= at <= 1.2:
                    biased.append(False)
                else:
                    biased.append(True)

    logger.debug("Finished bias assessment with result: {}".format(biased))

    # Split biased list in half to get the 5' and 3' ends
    biased_5end, biased_3end = biased[:int(len(biased)/2)],\
        biased[int(len(biased)/2):][::-1]

    logger.debug("Getting optimal trim range from biased list")
    trim_nt = [0, 0]
    # Assess number of nucleotides to clip at 5' end
    trim_nt[0] = get_trim_index(biased_5end)
    logger.debug("Optimal trim range at 5' end set to: {}".format(trim_nt[0]))
    # Assess number of nucleotides to clip at 3' end
    trim_nt[1] = len(biased) - get_trim_index(biased_3end)
    logger.debug("Optimal trim range at 3' end set to: {}".format(trim_nt[1]))

    return trim_nt


def get_sample_trim(p1_data, p2_data):
    """Get the optimal read trim range from data files of paired FastQ reads.

    Given the FastQC data report files for paired-end FastQ reads, this
    function will assess the optimal trim range for the 3' and 5' ends of
    the paired-end reads. This assessment will be based on the *'Per sequence
    GC content'*.

    Parameters
    ----------
    p1_data: str
        Path to FastQC data report file from pair 1
    p2_data: str
        Path to FastQC data report file from pair 2

    Returns
    -------
    optimal_5trim: int
        Optimal trim index for the 5' end of the reads
    optima_3trim: int
        Optimal trim index for the 3' end of the reads

    See Also
    --------
    trim_range

    """

    sample_ranges = [trim_range(x) for x in [p1_data, p2_data]]

    # Get the optimal trim position for 5' end
    optimal_5trim = max([x[0] for x in sample_ranges])
    # Get optimal trim position for 3' end
    optimal_3trim = min([x[1] for x in sample_ranges])

    return optimal_5trim, optimal_3trim


def get_summary(summary_file):
    """Parses a FastQC summary report file and returns it as a dictionary.

    This function parses a typical FastQC summary report file, retrieving
    only the information on the first two columns. For instance, a line could
    be::

        'PASS	Basic Statistics	SH10762A_1.fastq.gz'

    This parser will build a dictionary with the string in the second column
    as a key and the QC result as the value. In this case, the returned
    ``dict`` would be something like::

        {"Basic Statistics": "PASS"}

    Parameters
    ----------
    summary_file: str
        Path to FastQC summary report.

    Returns
    -------
    summary_info: :py:data:`OrderedDict`
        Returns the information of the FastQC summary report as an ordered
        dictionary, with the categories as strings and the QC result as values.

    """

    summary_info = OrderedDict()
    logger.debug("Retrieving summary information from file: {}".format(
        summary_file))

    with open(summary_file) as fh:
        for line in fh:
            # Skip empty lines
            if not line.strip():
                continue
            # Populate summary info
            fields = [x.strip() for x in line.split("\t")]
            summary_info[fields[1]] = fields[0]

    logger.debug("Retrieved summary information from file: {}".format(
        summary_info))

    return summary_info


def check_summary_health(summary_file, **kwargs):
    """Checks the health of a sample from the FastQC summary file.

    Parses the FastQC summary file and tests whether the sample is good
    or not. There are four categories that cannot fail, and two that
    must pass in order for the sample pass this check. If the sample fails
    the quality checks, a list with the failing categories is also returned.

    Categories that cannot fail::

        fail_sensitive = [
            "Per base sequence quality",
            "Overrepresented sequences",
            "Sequence Length Distribution",
            "Per sequence GC content"
        ]

    Categories that must pass::

        must_pass = [
            "Per base N content",
            "Adapter Content"
        ]

    Parameters
    ----------
    summary_file: str
        Path to FastQC summary file.

    Returns
    -------
    x : bool
        Returns ``True`` if the sample passes all tests. ``False`` if not.
    summary_info : list
        A list with the FastQC categories that failed the tests. Is empty
        if the sample passes all tests.
    """

    # Store the summary categories that cannot fail. If they fail, do not
    # proceed with this sample
    fail_sensitive = kwargs.get("fail_sensitive", [
        "Per base sequence quality",
        "Overrepresented sequences",
        "Sequence Length Distribution",
        "Per sequence GC content"
    ])
    logger.debug("Fail sensitive categories: {}".format(fail_sensitive))

    # Store summary categories that must pass. If they do not, do not proceed
    # with that sample
    must_pass = kwargs.get("must_pass", [
        "Per base N content",
        "Adapter Content"
    ])
    logger.debug("Must pass categories: {}".format(must_pass))

    warning_fail_sensitive = kwargs.get("warning_fail_sensitive", [
        "Per base sequence quality",
        "Overrepresented sequences",

    ])

    warning_must_pass = kwargs.get("warning_must_pass", [
        "Per base sequence content"
    ])

    # Get summary dictionary
    summary_info = get_summary(summary_file)

    # This flag will change to False if one of the tests fails
    health = True
    # List of failing categories
    failed = []
    # List of warning categories
    warning = []

    for cat, test in summary_info.items():

        logger.debug("Assessing category {} with result {}".format(cat, test))

        # FAILURES
        # Check for fail sensitive
        if cat in fail_sensitive and test == "FAIL":
            health = False
            failed.append("{}:{}".format(cat, test))
            logger.error("Category {} failed a fail sensitive "
                         "category".format(cat))

        # Check for must pass
        if cat in must_pass and test != "PASS":
            health = False
            failed.append("{}:{}".format(cat, test))
            logger.error("Category {} failed a must pass category".format(
                cat))

        # WARNINGS
        # Check for fail sensitive
        if cat in warning_fail_sensitive and test == "FAIL":
            warning.append("Failed category: {}".format(cat))
            logger.warning("Category {} flagged at a fail sensitive "
                           "category".format(cat))

        if cat in warning_must_pass and test != "PASS":
            warning.append("Did not pass category: {}".format(cat))
            logger.warning("Category {} flagged at a must pass "
                           "category".format(cat))

    # Passed all tests
    return health, failed, warning


@MainWrapper
def main(sample_id, result_p1, result_p2, opts):
    """Main executor of the fastqc_report template.

    If the "--ignore-tests" option is present in the ``opts`` argument,
    the health check of the sample will be bypassed, and it will pass the
    check. This option is used in the first run of FastQC. In the second
    run (after filtering with trimmomatic) this option is not provided and
    the samples are submitted to a health check before proceeding in the
    pipeline.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    result_p1 : list
        Two element list containing the path to the FastQC report files to
        the first FastQ pair.
        The first must be the nucleotide level report and the second the
        categorical report.
    result_p2: list
        Two element list containing the path to the FastQC report files to
        the second FastQ pair.
        The first must be the nucleotide level report and the second the
        categorical report.
    opts : list
        List of arbitrary options. See `Expected input`_.

    """

    logger.info("Starting fastqc report")
    json_dic = {}

    with open("{}_trim_report".format(sample_id), "w") as trep_fh, \
            open("optimal_trim", "w") as trim_fh, \
            open("{}_status_report".format(sample_id), "w") as rep_fh, \
            open(".status", "w") as status_fh, \
            open(".warning", "w") as warn_fh, \
            open(".fail", "w") as fail_fh, \
            open(".report.json", "w") as report_fh:

        # Perform health check according to the FastQC summary report for
        # each pair. If both pairs pass the check, send the 'pass' information
        # to the 'fastqc_health' channel. If at least one fails, send the
        # summary report.
        if "--ignore-tests" not in opts:

            # Get reports for each category in json format
            json_dic = write_json_report(sample_id, result_p1[0],
                                         result_p2[0])

            logger.info("Performing FastQ health check")
            for p, fastqc_summary in enumerate([result_p1[1], result_p2[1]]):

                logger.debug("Checking files: {}".format(fastqc_summary))
                # Get the boolean health variable and a list of failed
                # categories, if any
                health, f_cat, warnings = check_summary_health(fastqc_summary)
                logger.debug("Health checked: {}".format(health))
                logger.debug("Failed categories: {}".format(f_cat))

                # Write any warnings
                if warnings:
                    json_dic["warnings"] = [{
                        "sample": sample_id,
                        "table": "qc",
                        "value": []
                    }]
                    for w in warnings:
                        warn_fh.write("{}\\n".format(w))
                        json_dic["warnings"][0]["value"].append(w)

                # Rename category summary file to the channel that will publish
                # The results
                output_file = "{}_{}_summary.txt".format(sample_id, p)
                os.rename(fastqc_summary, output_file)
                logger.debug("Setting summary file name to {}".format(
                    output_file))

                # If one of the health flags returns False, send the summary
                # report through the status channel
                if not health:
                    fail_msg = "Sample failed quality control checks:" \
                               " {}".format(",".join(f_cat))
                    logger.warning(fail_msg)
                    fail_fh.write(fail_msg)
                    json_dic["fail"] = [{
                        "sample": sample_id,
                        "table": "qc",
                        "value": [fail_msg]
                    }]
                    report_fh.write(
                        json.dumps(json_dic, separators=(",", ":")))
                    status_fh.write("fail")
                    trim_fh.write("fail")
                    rep_fh.write("{}, {}\\n".format(sample_id, ",".join(f_cat)))
                    trep_fh.write("{},fail,fail\\n".format(sample_id))

                    return

            logger.info("Sample passed quality control checks")

        status_fh.write("pass")
        rep_fh.write("{}, pass\\n".format(sample_id))

        logger.info("Assessing optimal trim range for sample")
        # Get optimal trimming range for sample, based on the per base sequence
        # content
        optimal_trim = get_sample_trim(result_p1[0], result_p2[0])
        logger.info("Optimal trim range set to: {}".format(optimal_trim))
        trim_fh.write("{}".format(" ".join([str(x) for x in optimal_trim])))

        trep_fh.write("{},{},{}\\n".format(sample_id, optimal_trim[0],
                                           optimal_trim[1]))

        # The json dict report is only populated when the FastQC quality
        # checks are performed, that is, when the --ignore-tests option
        # is not provide
        if json_dic:
            report_fh.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == '__main__':

    main(SAMPLE_ID, RESULT_P1, RESULT_P2, OPTS)
