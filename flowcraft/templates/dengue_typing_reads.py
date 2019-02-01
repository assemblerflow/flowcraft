#!/usr/bin/env python3

"""
Purpose
-------

This module intends to type DENV genome assembly with seqTyping
(mapping mode)

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
    FASTQ_PAIR = '$fastq_pair'.split()
    REFERENCE = '$reference'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))
    logger.debug("REFERENCE: {}".format(REFERENCE))


def __get_version_seq_typing():

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
    for ch in ['/', '`', '*', '{', '}', '[', ']', '(', ')', '#', '+', '-', '.', '!', '\$', ':', '|']:
        text = text.replace(ch, "_")
    return text


def getSequence(ref, fasta):
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
    with open(file, "r") as typing_report:
        lines = typing_report.readlines()
    return lines[1].split('\\t')[3]


def getType(file):
    with open(file, "r") as result:
        return result.readline().strip()


def getConsesusSequence(best_reference, consensus, sample_id):
    gb_ID = best_reference.split('|')[0].replace(":", "_")
    fh_consensus = open(consensus, "r")

    entry = (x[1] for x in groupby(fh_consensus, lambda line: line[0] == ">"))

    for header in entry:

        headerStr = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in entry.__next__())

        if gb_ID in headerStr:
            with open(sample_id + '_consensus.fasta', "w") as output_file:
                output_file.write(
                    ">" + sample_id + "_consensus_" +
                    best_reference.split("_")[0] + "\\n" + seq.upper() + "\\n")

    fh_consensus.close()


def getScore(file):

    score = "fail"

    with open(file, "r") as typing_report:
        lines = typing_report.readlines()

        sequence_covered = lines[1].split("\\t")[4]
        sequence_identity = lines[1].split("\\t")[6]

    #TODO - score

    return score


@MainWrapper
def main(sample_id, assembly, fastq_pair, reference):
    """Main executor of the dengue_typing template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    assembly : str
        Assembly file.
    fastq_pair: list
        FastQ files
    reference: str
        String stating is the reference genome is to be recovered"""

    json_report = {}

    st_version = __get_version_seq_typing()

    cli = ["seq_typing.py",
           "reads",
           "--org", "Dengue", "virus",
           "-j", "${task.cpus}",
           "--doNotRemoveConsensus",
           "-f", fastq_pair[0], fastq_pair[1]]

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

        best_reference = get_reference_header("seq_typing.report_types.tab")

        if typing_result != "NT":
            logger.info("Getting consensus sequenceq")
            getConsesusSequence(best_reference,
                                "rematch/1_GenotypesDENV_14-05-18.headers_renamed.fasta_0/sample.noMatter.fasta",
                                sample_id)

            # TODO
            #confidence_score = getScore("seq_typing.report_types.tab")

        else:
            logger.error("Failed to obtain a close reference sequence in read mode. No consensus sequence is obtained.")
            with open(".status", "w") as status:
                status.write("fail")
            sys.exit(1)

    else:
        logger.error("Failed to run seq_typing for Dengue Virus.")
        with open(".status", "w") as status:
            status.write("fail")
        sys.exit(1)

    if reference == "true":
        reference_name = getSequence(best_reference,
                                     "/NGStools/seq_typing/seqtyping/reference_sequences/dengue_virus/1_GenotypesDENV_14-05-18.fasta")

        json_report = {'tableRow': [{
            'sample': sample_id,
            'data': [
                 {'header': 'seqtyping',
                  'value': typing_result,
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
                 'table': 'typing'}
            ]}],
            'metadata': [
                {'sample': sample_id,
                 'treeData': typing_result,
                 'column': 'typing'}]}

    # Add information to dotfiles
    with open(".report.json", "w") as report, \
            open(".status", "w") as status, \
            open(".version", "w") as version:
        report.write(json.dumps(json_report, separators=(",", ":")))
        status.write("pass")
        version.write(st_version)


if __name__ == '__main__':

    main(SAMPLE_ID, ASSEMBLY, FASTQ_PAIR, REFERENCE)
