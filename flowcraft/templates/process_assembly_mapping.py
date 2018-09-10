#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to process the coverage report from the
:py:class:`assembly_mapping` process.

TODO: Better purpose

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``assembly`` : Fasta assembly file.
    - e.g.: ``'SH10761A.assembly.fasta'``
- ``coverage`` : TSV file with the average coverage for each assembled contig.
    - e.g.: ``'coverage.tsv'``
- ``coverage_bp`` : TSV file with the coverage for each assembled bp.
    - e.g.: ``'coverage.tsv'``
- ``bam_file`` : BAM file with the alignment of reads to the genome.
    - e.g.: ``'sorted.bam'``
- ``opts`` : List of options for processing assembly mapping output.
    1. Minimum coverage for assembled contigs. Can be``auto``.
        - e.g.: ``'auto'`` or ``'10'``
    2. Maximum number of contigs.
        - e.g.: '100'
- ``gsize``: Expected genome size.
    - e.g.: ``'2.5'``

Generated output
----------------
- ``${sample_id}_filtered.assembly.fasta`` : Filtered assembly file in Fasta \
    format.
    - e.g.: ``'SampleA_filtered.assembly.fasta'``
- ``filtered.bam`` : BAM file with the same filtering as the assembly file.
    - e.g.: ``filtered.bam``


Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "09022018"
__template__ = "process_assembly_mapping-nf"

import os
import json
import shutil
import subprocess

from subprocess import PIPE
from collections import OrderedDict

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


def __get_version_samtools():

    try:
        cli = ["samtools", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout = p.communicate()[0]

        version = stdout.splitlines()[0].split()[1].decode("utf8")
    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return {
        "program": "Samtools",
        "version": version
    }


def __get_version_bowtie2():

    try:
        cli = ["bowtie2", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout = p.communicate()[0]

        version = stdout.splitlines()[0].split()[-1].decode("utf8")
    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return {
        "program": "Bowtie2",
        "version": version
    }


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLY_FILE = '$assembly'
    COVERAGE_FILE = '$coverage'
    COVERAGE_BP_FILE = '$coverage_bp'
    BAM_FILE = '$bam_file'
    OPTS = [x.strip() for x in '$opts'.strip("[]").split(",")]
    GSIZE = float('$gsize')
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLY_FILE: {}".format(ASSEMBLY_FILE))
    logger.debug("COVERAGE_FILE: {}".format(COVERAGE_FILE))
    logger.debug("COVERAGE_BP_FILE: {}".format(COVERAGE_BP_FILE))
    logger.debug("BAM_FILE: {}".format(BAM_FILE))
    logger.debug("MIN_ASSEMBLY_COVERAGE: {}".format(OPTS))
    logger.debug("GSIZE: {}".format(GSIZE))


def parse_coverage_table(coverage_file):
    """Parses a file with coverage information into objects.

    This function parses a TSV file containing coverage results for
    all contigs in a given assembly and will build an ``OrderedDict``
    with the information about their coverage and length.  The length
    information is actually gathered from the contig header using a
    regular expression that assumes the usual header produced by Spades::

        contig_len = int(re.search("length_(.+?)_", line).group(1))

    Parameters
    ----------
    coverage_file : str
        Path to TSV file containing the coverage results.

    Returns
    -------
    coverage_dict : OrderedDict
        Contains the coverage and length information for each contig.
    total_size : int
        Total size of the assembly in base pairs.
    total_cov : int
        Sum of coverage values across all contigs.
    """

    # Stores the correspondence between a contig and the corresponding coverage
    # e.g.: {"contig_1": {"cov": 424} }
    coverage_dict = OrderedDict()
    # Stores the total coverage
    total_cov = 0

    with open(coverage_file) as fh:
        for line in fh:
            # Get contig and coverage
            contig, cov = line.strip().split()
            coverage_dict[contig] = {"cov": int(cov)}
            # Add total coverage
            total_cov += int(cov)
            logger.debug("Processing contig '{}' with coverage '{}'"
                         "".format(contig, cov))

    return coverage_dict, total_cov


def filter_assembly(assembly_file, minimum_coverage, coverage_info,
                    output_file):
    """Generates a filtered assembly file.

    This function generates a filtered assembly file based on an original
    assembly and a minimum coverage threshold.

    Parameters
    ----------
    assembly_file : str
        Path to original assembly file.
    minimum_coverage : int or float
        Minimum coverage required for a contig to pass the filter.
    coverage_info : OrderedDict or dict
        Dictionary containing the coverage information for each contig.
    output_file : str
        Path where the filtered assembly file will be generated.

    """

    # This flag will determine whether sequence data should be written or
    # ignored because the current contig did not pass the minimum
    # coverage threshold
    write_flag = False

    with open(assembly_file) as fh, open(output_file, "w") as out_fh:

        for line in fh:
            if line.startswith(">"):
                # Reset write_flag
                write_flag = False
                # Get header of contig
                header = line.strip()[1:]
                # Check coverage for current contig
                contig_cov = coverage_info[header]["cov"]
                # If the contig coverage is above the threshold, write to
                # output filtered assembly
                if contig_cov >= minimum_coverage:
                    write_flag = True
                    out_fh.write(line)

            elif write_flag:
                out_fh.write(line)


def filter_bam(coverage_info, bam_file, min_coverage, output_bam):
    """Uses Samtools to filter a BAM file according to minimum coverage

    Provided with a minimum coverage value, this function will use Samtools
    to filter a BAM file. This is performed to apply the same filter to
    the BAM file as the one applied to the assembly file in
    :py:func:`filter_assembly`.

    Parameters
    ----------
    coverage_info : OrderedDict or dict
        Dictionary containing the coverage information for each contig.
    bam_file : str
        Path to the BAM file.
    min_coverage : int
        Minimum coverage required for a contig to pass the filter.
    output_bam : str
        Path to the generated filtered BAM file.
    """

    # Get list of contigs that will be kept
    contig_list = [x for x, vals in coverage_info.items()
                   if vals["cov"] >= min_coverage]

    cli = [
        "samtools",
        "view",
        "-bh",
        "-F",
        "4",
        "-o",
        output_bam,
        "-@",
        "1",
        bam_file,
    ]

    cli += contig_list

    logger.debug("Runnig samtools view subprocess with command: {}".format(
        cli))

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # Attempt to decode STDERR output from bytes. If unsuccessful, coerce to
    # string
    try:
        stderr = stderr.decode("utf8")
        stdout = stdout.decode("utf8")
    except (UnicodeDecodeError, AttributeError):
        stderr = str(stderr)
        stdout = str(stdout)

    logger.info("Finished samtools view subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished samtools view subprocesswith STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished samtools view with return code: {}".format(
        p.returncode))

    if not p.returncode:
        # Create index
        cli = [
            "samtools",
            "index",
            output_bam
        ]

        logger.debug("Runnig samtools index subprocess with command: "
                     "{}".format(cli))

        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        try:
            stderr = stderr.decode("utf8")
            stdout = stdout.decode("utf8")
        except (UnicodeDecodeError, AttributeError):
            stderr = str(stderr)
            stdout = str(stdout)

        logger.info("Finished samtools index subprocess with STDOUT:\\n"
                    "======================================\\n{}".format(
            stdout))
        logger.info("Fished samtools index subprocesswith STDERR:\\n"
                    "======================================\\n{}".format(
            stderr))
        logger.info("Finished samtools index with return code: {}".format(
            p.returncode))


def check_filtered_assembly(coverage_info, coverage_bp, minimum_coverage,
                            genome_size, contig_size, max_contigs,
                            sample_id):
    """Checks whether a filtered assembly passes a size threshold

    Given a minimum coverage threshold, this function evaluates whether an
    assembly will pass the minimum threshold of ``genome_size * 1e6 * 0.8``,
    which means 80% of the expected genome size or the maximum threshold
    of ``genome_size * 1e6 * 1.5``, which means 150% of the expected genome
    size. It will issue a warning if any of these thresholds is crossed.
    In the case of an expected genome size below 80% it will return False.

    Parameters
    ----------
    coverage_info : OrderedDict or dict
        Dictionary containing the coverage information for each contig.
    coverage_bp : dict
        Dictionary containing the per base coverage information for each
        contig. Used to determine the total number of base pairs in the
        final assembly.
    minimum_coverage : int
        Minimum coverage required for a contig to pass the filter.
    genome_size : int
        Expected genome size.
    contig_size : dict
        Dictionary with the len of each contig. Contig headers as keys and
        the corresponding lenght as values.
    max_contigs : int
        Maximum threshold for contig number. A warning is issued if this
        threshold is crossed.
    sample_id : str
        Id or name of the current sample

    Returns
    -------
    x : bool
        True if the filtered assembly size is higher than 80% of the
        expected genome size.

    """

    # Get size of assembly after filtering contigs below minimum_coverage
    assembly_len = sum([v for k, v in contig_size.items()
                        if coverage_info[k]["cov"] >= minimum_coverage])
    logger.debug("Assembly length after filtering for minimum coverage of"
                 " {}: {}".format(minimum_coverage, assembly_len))
    # Get number of contigs after filtering
    ncontigs = len([x for x in coverage_info.values()
                    if x["cov"] >= minimum_coverage])
    logger.debug("Number of contigs: {}".format(ncontigs))
    # Get number of bp after filtering
    filtered_contigs = [k for k, v in coverage_info.items()
                        if v["cov"] >= minimum_coverage]
    logger.debug("Filtered contigs for minimum coverage of "
                 "{}: {}".format(minimum_coverage, filtered_contigs))
    total_assembled_bp = sum([sum(coverage_bp[x]) for x in filtered_contigs
                              if x in coverage_bp])
    logger.debug("Total number of assembled base pairs:"
                 "{}".format(total_assembled_bp))

    warnings = []
    fails = []
    health = True

    with open(".warnings", "w") as warn_fh, \
            open(".report.json", "w") as json_report:

        logger.debug("Checking assembly size after filtering : {}".format(
            assembly_len))

        # If the filtered assembly size is above the 150% genome size
        # threshold, issue a warning
        if assembly_len > genome_size * 1e6 * 1.5:
            warn_msg = "Assembly size ({}) smaller than the maximum" \
                       " threshold of 150% of expected genome size.".format(
                            assembly_len)
            logger.warning(warn_msg)
            warn_fh.write(warn_msg)
            fails.append("Large_genome_size_({})".format(assembly_len))

        # If the number of contigs in the filtered assembly size crosses the
        # max_contigs threshold, issue a warning
        logger.debug("Checking number of contigs: {}".format(
                len(coverage_info)))
        contig_threshold = max_contigs * genome_size / 1.5
        if ncontigs > contig_threshold:
            warn_msg = "The number of contigs ({}) exceeds the threshold of " \
                       "100 contigs per 1.5Mb ({})".format(
                            ncontigs, round(contig_threshold, 1))
            logger.warning(warn_msg)
            warn_fh.write(warn_msg)
            warnings.append(warn_msg)

        # If the filtered assembly size falls below the 80% genome size
        # threshold, fail this check and return False
        if assembly_len < genome_size * 1e6 * 0.8:
            warn_msg = "Assembly size smaller than the minimum" \
                       " threshold of 80% of expected genome size: {}".format(
                            assembly_len)
            logger.warning(warn_msg)
            warn_fh.write(warn_msg)
            fails.append("Small_genome_size_({})".format(assembly_len))
            assembly_len = sum([v for v in contig_size.values()])
            total_assembled_bp = sum(
                [sum(coverage_bp[x]) for x in coverage_info if x in
                 coverage_bp])
            logger.debug("Assembly length without coverage filtering: "
                         "{}".format(assembly_len))
            logger.debug("Total number of assembled base pairs without"
                         " filtering: {}".format(total_assembled_bp))

            health = False

        json_dic = {
            "plotData": [{
                "sample": sample_id,
                "data": {
                    "sparkline": total_assembled_bp
                }
            }]
        }

        if warnings:
            json_dic["warnings"] = [{
                "sample": sample_id,
                "table": "assembly",
                "value": warnings
            }]
        if fails:
            json_dic["fail"] = [{
                "sample": sample_id,
                "table": "assembly",
                "value": [fails]
            }]

        json_report.write(json.dumps(json_dic, separators=(",", ":")))

    return health


def get_coverage_from_file(coverage_file):
    """

    Parameters
    ----------
    coverage_file

    Returns
    -------

    """

    contig_coverage = {}

    with open(coverage_file) as fh:
        for line in fh:

            fields = line.strip().split()

            # Get header
            header = fields[0]
            coverage = int(fields[2])

            if header not in contig_coverage:
                contig_coverage[header] = [coverage]
            else:
                contig_coverage[header].append(coverage)

    return contig_coverage


def evaluate_min_coverage(coverage_opt, assembly_coverage, assembly_size):
    """ Evaluates the minimum coverage threshold from the value provided in
    the coverage_opt.

    Parameters
    ----------
    coverage_opt : str or int or float
        If set to "auto" it will try to automatically determine the coverage
        to 1/3 of the assembly size, to a minimum value of 10. If it set
        to a int or float, the specified value will be used.
    assembly_coverage : int or float
        The average assembly coverage for a genome assembly. This value
        is retrieved by the `:py:func:parse_coverage_table` function.
    assembly_size : int
        The size of the genome assembly. This value is retrieved by the
        `py:func:get_assembly_size` function.

    Returns
    -------
    x: int
        Minimum coverage threshold.

    """

    if coverage_opt == "auto":
        # Get the 1/3 value of the current assembly coverage
        min_coverage = (assembly_coverage / assembly_size) * .3
        logger.info("Minimum assembly coverage automatically set to: "
                    "{}".format(min_coverage))
        # If the 1/3 coverage is lower than 10, change it to the minimum of
        # 10
        if min_coverage < 10:
            logger.info("Minimum assembly coverage cannot be set to lower"
                        " that 10. Setting to 10")
            min_coverage = 10
    else:
        min_coverage = int(coverage_opt)
        logger.info("Minimum assembly coverage manually set to: {}".format(
            min_coverage))

    return min_coverage


def get_assembly_size(assembly_file):
    """Returns the number of nucleotides and the size per contig for the
    provided assembly file path

    Parameters
    ----------
    assembly_file : str
        Path to assembly file.

    Returns
    -------
    assembly_size : int
        Size of the assembly in nucleotides
    contig_size : dict
        Length of each contig (contig name as key and length as value)

    """

    assembly_size = 0
    contig_size = {}
    header = ""

    with open(assembly_file) as fh:
        for line in fh:

            # Skip empty lines
            if line.strip() == "":
                continue

            if line.startswith(">"):
                header = line.strip()[1:]
                contig_size[header] = 0

            else:
                line_len = len(line.strip())
                assembly_size += line_len
                contig_size[header] += line_len

    return assembly_size, contig_size


@MainWrapper
def main(sample_id, assembly_file, coverage_file, coverage_bp_file, bam_file,
         opts, gsize):
    """Main executor of the process_assembly_mapping template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    assembly_file : str
        Path to assembly file in Fasta format.
    coverage_file : str
        Path to TSV file with coverage information for each contig.
    coverage_bp_file : str
        Path to TSV file with coverage information for each base.
    bam_file : str
        Path to BAM file.
    opts : list
        List of options for processing assembly mapping.
    gsize : int
        Expected genome size

    """

    min_assembly_coverage, max_contigs = opts

    logger.info("Starting assembly mapping processing")

    # Get coverage info, total size and total coverage from the assembly
    logger.info("Parsing coverage table")
    coverage_info, a_cov = parse_coverage_table(coverage_file)
    a_size, contig_size = get_assembly_size(assembly_file)
    logger.info("Assembly processed with a total size of '{}' and coverage"
                " of '{}'".format(a_size, a_cov))
    # Get number of assembled bp after filters
    logger.info("Parsing coverage per bp table")
    coverage_bp_data = get_coverage_from_file(coverage_bp_file)

    # Assess the minimum assembly coverage
    min_coverage = evaluate_min_coverage(min_assembly_coverage, a_cov, a_size)

    # Check if filtering the assembly using the provided min_coverage will
    # reduce the final bp number to less than 80% of the estimated genome
    # size.
    # If the check below passes with True, then the filtered assembly
    # is above the 80% genome size threshold.
    filtered_assembly = "{}_filt.fasta".format(
        os.path.splitext(assembly_file)[0])
    filtered_bam = "filtered.bam"
    logger.info("Checking filtered assembly")
    if check_filtered_assembly(coverage_info, coverage_bp_data, min_coverage,
                               gsize, contig_size, int(max_contigs),
                               sample_id):
        # Filter assembly contigs based on the minimum coverage.
        logger.info("Filtered assembly passed minimum size threshold")
        logger.info("Writting filtered assembly")
        filter_assembly(assembly_file, min_coverage, coverage_info,
                        filtered_assembly)
        logger.info("Filtering BAM file according to saved contigs")
        filter_bam(coverage_info, bam_file, min_coverage, filtered_bam)
    # Could not filter the assembly as it would drop below acceptable
    # length levels. Copy the original assembly to the output assembly file
    # for compliance with the output channel
    else:
        shutil.copy(assembly_file, filtered_assembly)
        shutil.copy(bam_file, filtered_bam)
        shutil.copy(bam_file + ".bai", filtered_bam + ".bai")

    with open(".status", "w") as status_fh:
        status_fh.write("pass")


if __name__ == '__main__':

    main(SAMPLE_ID, ASSEMBLY_FILE, COVERAGE_FILE, COVERAGE_BP_FILE,
         BAM_FILE, OPTS, GSIZE)
