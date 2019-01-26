#!/usr/bin/env python3

"""
Purpose
-------

This module is intended execute mlst on Fasta files.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fasta_file`` : Fasta file paths.
    - e.g.: ``'SampleA.fasta'``
- ``mlstSpecies`` : Expected species

Generated output
----------------


Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "09012019"
__template__ = "mlst-nf"

import os
import sys
import subprocess
import json

from itertools import groupby
from subprocess import PIPE

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLY = '$assembly'
    EXPECTED_SPECIES = '$expected_species'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("EXPECTED_SPECIES: {}".format(EXPECTED_SPECIES))


def chunkstring(string, length):
    """
    Divides sequences in a multifasta file to the provided length.

    Parameters
    ----------
    string: Str
        Sequence to divide.
    length: int
        maximum sequence length.
    """
    return (string[0 + i:length + i] for i in range(0, len(string), length))


def get_species_scheme_map_version(mlst_folder):
    """
    Since release v2.16.1, the file that maps the schema genus
    to the species name changed from "species_scheme_map" to
    "scheme_soecies_map.tab". This method determines which version
    is in the container provided and returns it. If no file is found,
    it terminates the code execution.

    Parameters
    ----------
    mlst_folder: str
        Path to the mlst source code

    Returns
    -------
    mlst_db_path: str
        Path to the mlst database containing the species and respective schema
    species_scheme_map_version: int
        Version on the mlst scheme database
    """

    mlst_db_path_version1 = os.path.join(os.path.dirname(os.path.dirname(mlst_folder)), 'db', 'species_scheme_map.tab')

    mlst_db_path_version2 = os.path.join(os.path.dirname(os.path.dirname(mlst_folder)), 'db', 'scheme_species_map.tab')

    if os.path.isfile(mlst_db_path_version1):
        return mlst_db_path_version1, 1

    elif os.path.isfile(mlst_db_path_version2):
        return mlst_db_path_version2, 2

    else:
        logger.error("Species_scheme_map not found.")
        sys.exit(1)


def set_species_scheme_map_variables(list_values, species_scheme_map_version):
    """
    Depending on the version of the mlst database containing
    the species and respective schema, it retrieves the entries for the
    genus, species and scheme name.

    Parameters
    ----------
    list_values: list
        line, in list form, of the mlst scheme database
    species_scheme_map_version: int
        Version of the mlst dabase.

    Returns
    -------
    val_genus: str
        genus in line
    val_species: str
        species in line
    Val_scheme: str
        scheme name in line
    """

    if species_scheme_map_version == 1:
        val_genus = list_values[0]
        val_species = list_values[1]
        val_scheme = list_values[2]
    elif species_scheme_map_version == 2:
        val_genus = list_values[1]
        val_species = list_values[2]
        val_scheme = list_values[0]

    return val_genus, val_species, val_scheme


def parse_species_scheme_map(species_splited, mlst_db_path, species_scheme_map_version):
    """
    Parses the mlst scheme database and returns the full scheme for the species
    and the genus.

    Parameters
    ----------
    species_splited: list
        List with the genus and specific epithet for the expected species
    mlst_db_path: str
        Path to the mlst database containing the species and respective schema
    species_scheme_map_version: int
        Version on the mlst database

    Returns
    -------
    scheme: str
        mlst scheme name
    genus_mlst_scheme: str
        genus name for the mlst scheme
    """
    scheme = 'unknown'
    genus_mlst_scheme = None

    with open(mlst_db_path, 'rtU') as reader:
        for line in reader:
            scheme_line = line.splitlines()[0]

            if scheme_line and not scheme_line.startswith('#'):
                scheme_line = scheme_line.lower().split('\\t')
                scheme_line_data = [scheme_line[i].split(' ')[0] for i in range(0, len(scheme_line))]
                val_genus, val_species, val_scheme = set_species_scheme_map_variables(scheme_line_data,
                                                                                      species_scheme_map_version)
                # checking if genus from expected species and scheme match
                if val_genus == species_splited[0]:

                    # if the scheme is not species specific (works for genus), the genus is set as the scheme name
                    if val_species == '':
                        genus_mlst_scheme = val_scheme

                    elif val_species == species_splited[1]:
                        scheme = val_scheme

        if scheme == 'unknown' and genus_mlst_scheme is not None:
            scheme = genus_mlst_scheme

    return scheme, genus_mlst_scheme


def clean_novel_alleles(novel_alleles, scheme_mlst, profile):
    """
    Clean the fasta file with the novel alleles produced by mlst

    Parameters
    ----------
    novel_alleles : str
        Path for fasta file containing the novel alleles
    scheme_mlst : str
        MLST schema found by mlst
    profile : list
        List of strings with the profile found

    Returns
    -------
    """
    unknown_genes = []

    #get novel alleles
    for gene_allele in profile:
        print(gene_allele)
        gene = gene_allele.split('(')[0]
        try:
            allele = gene_allele.split('(')[1].rstrip(')')
            if allele.startswith('~'):
                unknown_genes.append(gene)
        except IndexError as e:
            logger.warning("WARNING: Found unexpected formatting on mlst profile {}".format(e))

    novel_alleles_keep = {}

    if unknown_genes:

        try:
            reader = open(novel_alleles, mode='rt', newline=None)

            fasta_iter = (g for k, g in groupby(reader, lambda x: x.startswith('>')))

            for header in fasta_iter:

                header = header.__next__()[1:].rstrip('\\r\\n')
                seq = ''.join(s.rstrip('\\r\\n') for s in fasta_iter.__next__())

                if header.startswith(scheme_mlst):
                    gene = header.split('.')[1].split('~')[0]

                    if gene in unknown_genes:
                        novel_alleles_keep[header] = seq
            reader.close()

            os.remove(novel_alleles)

            if novel_alleles_keep:
                with open(novel_alleles, 'wt') as writer:
                    for header, seq in novel_alleles_keep.items():
                        writer.write('>{}\\n'.format(header))
                        writer.write('\\n'.join(chunkstring(seq, 80)) + '\\n')
        except FileNotFoundError as e:
            logger.info("An unknown ST was found but no novel alleles fasta file was "
                        "produced by mlst software: {}".format(e))


def getScheme(species):
    """
    Get mlst scheme for the expected species.

    Parameters
    ----------
    species: str
        Expected species

    Returns
    -------
    scheme : str
        mlst scheme name
    species_genus : str
        genus of the expected species
    genus_mlst_scheme : str
        genus for the mlst scheme
    """
    cli = ['which', 'mlst']

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

    mlst_folder = os.path.abspath(os.path.realpath(stdout.splitlines()[0]))

    mlst_db_path, species_scheme_map_new = get_species_scheme_map_version(mlst_folder)

    scheme, genus_mlst_scheme = parse_species_scheme_map(species.lower().split(' '), mlst_db_path,
                                                         species_scheme_map_new)

    logger.info('MLST scheme found for {species}: {scheme}'.format(species=species, scheme=scheme))

    species_genus = species.lower().split(' ')[0]

    return scheme, species_genus, genus_mlst_scheme


def parse_stdout(stdout):
    """
    Parses mlst stdout to retrieve the hit's mlst
    scheme, st and profile

    Parameters
    ----------
    stdout: str
        mlst stdout

    Returns
    -------
    scheme_mlst : str
        mlst scheme name
    st : str
        mlst st for the sample
    profile : list
        list of strings containing the profile

    """

    mlst_data = stdout.splitlines()[0].split("\\t")

    scheme_mlst = mlst_data[1].split("_")[0]
    st = mlst_data[2]
    profile = mlst_data[3:]

    return scheme_mlst, st, profile


@MainWrapper
def main(sample_id, assembly, expected_species):
    """
    Main executor of the mlst template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    assembly : str
        Fasta file.
    expected_species : str
        Expected species

    """

    pass_qc = False

    novel_alleles = "{}_mlst_novel_alleles.fasta".format(sample_id)

    cli = [
        "mlst",
        "--novel",
        novel_alleles,
        assembly
    ]

    logger.debug("Running mlst subprocess with command: {}".format(cli))

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

    logger.info("Finished mlst subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished mlst subprocess with STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished mlst with return code: {}".format(
        p.returncode))

    if p.returncode == 0:

        with open("{}.mlst.txt".format(sample_id), "wt") as writer:
            writer.write(stdout)

        scheme_mlst, st, profile = parse_stdout(stdout)

        # In case it's an unkown ST, cleans the novel alleles file.
        if st == "-":
            clean_novel_alleles(novel_alleles=novel_alleles, scheme_mlst=scheme_mlst, profile=profile)
        else:
            # in case mlst fails to create the novel alleles file
            if os.path.isfile(novel_alleles):
                os.remove(novel_alleles)

        # if the expected_species is set to PASS, it bypasses species verification
        if expected_species == "PASS":
            pass_qc = True

        else:
            # returns the schema for the expected species, the genus for that species and the genus for the mlst schema
            expected_scheme, expected_species_genus, expected_mlst_scheme_genus = getScheme(expected_species)

            if scheme_mlst.split('_', 1)[0] == expected_scheme.split("_", 1)[0]:
                pass_qc = True

            # If the scheme is not the same as the expected species,
            else:
                if expected_scheme == "unknown":
                    pass_qc = True
                    if scheme_mlst != "-":
                        logger.warning("Found {} scheme for expected species".format(scheme_mlst))

                # in case of yersinia, it passes QC if the genus match as there's a scheme just for the genus,
                # and one for yersinia pseudotuberculosis
                elif expected_species_genus == 'yersinia' and expected_mlst_scheme_genus == 'yersinia':
                    pass_qc = True
                    logger.warning("Found a Yersinia scheme ({}), but it is different from what it was"
                                   " expected ({})".format(scheme_mlst, expected_scheme))
                else:
                    if expected_mlst_scheme_genus is not None and \
                            scheme_mlst == expected_scheme == expected_mlst_scheme_genus:
                        pass_qc = True
                    else:
                        logger.error("MLST scheme found ({}) and provided ({}) are not the same"
                                     .format(scheme_mlst, expected_scheme))

        # writing .report.json
        report_dic = {
            "expectedSpecies": expected_species,
            "species": scheme_mlst,
            "st": st,
            "tableRow": [
                {"sample": sample_id,
                 "data": [
                     {'header': 'MLST species',
                      'value': scheme_mlst,
                      'table': 'typing'
                      },
                     {'header': 'MLST ST',
                      'value': st,
                      'table': 'typing',
                      "columnBar": False
                      }
                 ]}
            ]
        }

        with open(".report.json", "w") as report_fh, \
                open(".status", "w") as status_fh:

            report_fh.write(json.dumps(report_dic, separators=(",", ":")))

            if pass_qc:
                status_fh.write("pass")
            else:
                status_fh.write("fail")
    else:
        logger.error("Sample {} did not run successfully".format(sample_id))
        with open(".status", "w") as status_fh:
            status_fh.write("error")


if __name__ == '__main__':

    main(SAMPLE_ID, ASSEMBLY, EXPECTED_SPECIES)
