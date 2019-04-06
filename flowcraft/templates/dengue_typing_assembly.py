#!/usr/bin/env python3

"""
Purpose
-------

This module intends to type DENV genome assembly with seqTyping (BLAST mode)

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fasta`` : A fasta file path.
    - e.g.: ``'SampleA.fasta'``
- ``fastq_pair`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'``

Generated output
----------------

-  The sample fasta file path or, if a complete ORF isn't obtained, a consesus sequence
-  The closest reference fasta file path
"""

__version__ = "0.0.2"
__build__ = "01022019"
__template__ = "dengue_typing-nf"

import json
import os
import sys
import subprocess
from subprocess import PIPE
from itertools import groupby
from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLY = '$assembly'
    REFERENCE = '$reference'
    RESULT = '$get_reference'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("REFERENCE: {}".format(REFERENCE))
    logger.debug("RESULT: {}".format(RESULT))


def __get_version_seq_typing():
    """
    Gets Seq_typing software version
    Returns
    -------
    version : str
        Seqtyping version"""

    try:
        cli = ["seq_typing.py", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout = p.communicate()[0]

        version = stdout.splitlines()[0].split()[-1].decode("utf8")
    except Exception as e:
        logger.debug(e)
        version = "undefined"

    return version

def replace_char(text):
    """
    Cleans the string from problematic chars

    Parameters
    ----------
    text : str
        String to clean"""

    for ch in ['/', '`', '*', '{', '}', '[', ']', '(', ')', '#', '+', '-', '.', '!', '\$', ':', '|']:
        text = text.replace(ch, "_")
    return text


def getSequence(ref, fasta):
    """
     Gets the fasta sequence from the Database with the header "ref"

     Parameters
     ----------
     ref : str
         Reference whose sequence needs to be fetched
     fasta: str
        Path to the multifasta"""

    fasta_header = ""

    fh_fasta = open(fasta, "r")
    entry = (x[1] for x in groupby(fh_fasta, lambda line: line[0] == ">"))

    for header in entry:
        headerStr = header.__next__()[1:].strip()

        seq = "".join(s.strip() for s in entry.__next__())

        if ref == headerStr.replace('>',''):
            filename = os.path.join(os.getcwd(), ref.replace('/','_').split('|')[0])
            fasta_header = replace_char(headerStr)

            with open(filename + '.fa', "w") as output_file:
                output_file.write(">" + fasta_header + "\\n" + seq.upper() + "\\n")

    fh_fasta.close()
    return fasta_header


def get_reference_header(file):
    """
    Gets the header for the closest reference from the seqtyping report

    Parameters
    ----------
    file: str
     Path to the seqtyping report"""

    with open(file, "r") as typing_report:
        lines = typing_report.readlines()
    return lines[1].split('\\t')[3]


def getType(file):
    """
    Gets the typing result from the seqtyping report

    Parameters
    ----------
    file: str
     Path to the seqtyping report"""

    with open(file, "r") as result:
        return result.readline().strip()


def getScore(file):
    """
    Method to write QC warnings based on the mapping statistics
    (sequence covered and identity)

    Parameters
    ----------
    file: str
     Path to the seqtyping report"""

    with open(file, "r") as typing_report:
        lines = typing_report.readlines()

        sequence_covered = float(lines[1].split("\\t")[4])
        sequence_identity = float(lines[1].split("\\t")[6])

        if sequence_covered < 70:
            logger.fail("Sequence coverage below 70% on the best hit.")
            with open(".fails", "w") as fails:
                fails.write("Sequence coverage below 70% on the best hit.")

        elif 90 > sequence_covered < 70:
            logger.warning("Sequence coverage lower than 90% on the best hit.")
            with open(".warnings", "w") as fails:
                fails.write("Sequence coverage below 70% on the best hit.")

        return sequence_identity, sequence_covered

@MainWrapper
def main(sample_id, assembly, reference, result):
    """Main executor of the dengue_typing template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    assembly : str
        Assembly file.
    fastq_pair: list
        FastQ files
    result: str
        String stating is the reference genome is to be recovered"""

    json_report = {}

    st_version = __get_version_seq_typing()

    cli = ["seq_typing.py",
           "assembly",
           "-b", os.path.join(os.getcwd(), reference),
           "-j", "${task.cpus}",
           "-f", assembly,
           "-t", "nucl"]

    logger.info("Runnig seq_typing subprocess with command: {}".format(cli))

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    try:
        stderr = stderr.decode("utf8")
        stdout = stdout.decode("utf8")
    except (UnicodeDecodeError, AttributeError):
        stderr = str(stderr)
        stdout = str(stdout)

    logger.info("Finished seq_typing index subprocess with STDOUT:\\n"
                "======================================\\n{}".format(
        stdout))
    logger.info("Fished seq_typing index subprocesswith STDERR:\\n"
                "======================================\\n{}".format(
        stderr))
    logger.info("Finished seq_typing index with return code: {}".format(
        p.returncode))

    if p.returncode == 0:

        typing_result = getType("seq_typing.report.txt")

        logger.info("Type found: {}".format(typing_result))

        # write appropriate QC dot files based on blast statistics
        identity = 0
        coverage = 0

        if typing_result != "NT":
            # write appropriate QC dot files based on blast statistics
            identity, coverage = getScore("seq_typing.report_types.tab")

            best_reference = get_reference_header("seq_typing.report_types.tab")

            reference_name = getSequence(best_reference, os.path.join(os.getcwd(), reference))

        else:
            logger.info("No typing information was obtained.")

        if result == "true":

            json_report = {'tableRow': [{
                'sample': sample_id,
                'data': [
                    {'header': 'seqtyping',
                     'value': typing_result,
                     'table': 'typing'},
                    {'header': 'Identity',
                     'value': round(identity,2),
                     'table': 'typing'},
                    {'header': 'Coverage',
                     'value': round(coverage,2),
                     'table': 'typing'},
                    {'header': 'Reference',
                     'value': reference_name.replace("gb_", "gb:").split("_")[0],
                     'table': 'typing'}
                ]}],
                'metadata': [
                    {'sample': sample_id,
                     'treeData': typing_result,
                     'column': 'typing'},
                    {'sample': reference_name,
                     'treeData': typing_result,
                     'column': 'typing'}]}

        else:
            json_report = {'tableRow': [{
                'sample': sample_id,
                'data': [
                    {'header': 'seqtyping',
                     'value': typing_result,
                     'table': 'typing'},
                    {'header': 'Identity',
                     'value': round(identity),
                     'table': 'typing'},
                    {'header': 'Coverage',
                     'value': round(coverage),
                     'table': 'typing'},
                    {'header': 'Reference',
                     'value': reference_name.replace("gb_", "gb:").split("_")[1],
                     'table': 'typing'}
                ]}],
                'metadata': [
                    {'sample': sample_id,
                     'treeData': typing_result,
                     'column': 'typing'}]}

    else:
        logger.error("Failed to run seq_typing for Dengue Virus.")
        with open(".status", "w") as status:
            status.write("fail")
        sys.exit(1)

    # Add information to dotfiles
    with open(".report.json", "w") as report, \
            open(".status", "w") as status, \
            open(".version", "w") as version:
        report.write(json.dumps(json_report, separators=(",", ":")))
        status.write("pass")
        version.write(st_version)


if __name__ == '__main__':

    main(SAMPLE_ID, ASSEMBLY, REFERENCE, RESULT)
