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
    is in the container provided and returns it.

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
    species_scheme_map_version = 1

    mlst_db_path = os.path.join(os.path.dirname(os.path.dirname(mlst_folder)), 'db', 'species_scheme_map.tab')

    if not os.path.isfile(mlst_db_path):
        mlst_db_path = os.path.join(os.path.dirname(os.path.dirname(mlst_folder)), 'db', 'scheme_species_map.tab')

        if not os.path.isfile(mlst_db_path):
            logger.error("Species_scheme_map not found.")
            sys.exit(1)

        else:
            species_scheme_map_version = 2

    return mlst_db_path, species_scheme_map_version


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
            line = line.splitlines()[0]

            if len(line) > 0:

                if not line.startswith('#'):
                    line = line.lower().split('\\t')
                    line = [line[i].split(' ')[0] for i in range(0, len(line))]
                    val_genus, val_species, val_scheme = set_species_scheme_map_variables(line, species_scheme_map_version)

                    if val_genus == species_splited[0]:

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
    for gene_allele in profile:
        gene = gene_allele.split('(')[0]
        try:
            allele = gene_allele.split('(')[1].rstrip(')')
            if allele.startswith('~'):
                unknown_genes.append(gene)
        except IndexError as e:
            print('WARNING: {}'.format(e))

    try:
        novel_alleles_keep = {}
        if len(unknown_genes) > 0:
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

        if len(novel_alleles_keep) > 0:
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

    novel_alleles = '{}_mlst_novel_alleles.fasta'.format(sample_id)

    cli = [
        'mlst',
        '--novel',
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

    with open(".status", "w") as fh:
        if p.returncode != 0:
            fh.write("error")
            return
        else:
            fh.write("pass")

    if p.returncode == 0:

        with open('{}.mlst.txt'.format(sample_id), 'wt') as writer:
            writer.write(stdout)

        # str
        scheme_mlst = stdout.splitlines()[0].split('\\t')[1].split('_')[0]
        # str
        st = stdout.splitlines()[0].split('\\t')[2]
        # list
        profile = stdout.splitlines()[0].split('\\t')[3:]

        # In case it's an unkown ST, cleans the novel alleles file.
        if st == '-':
            clean_novel_alleles(novel_alleles=novel_alleles, scheme_mlst=scheme_mlst, profile=profile)
        else:
            # in case mlst fails to create the novel alleles file
            if os.path.isfile(novel_alleles):
                os.remove(novel_alleles)

        if not expected_species == "PASS":

            # returns the schema for the expected species, the genus for that species and the genus for the mlst schema
            scheme, species_genus, mlst_scheme_genus = getScheme(expected_species)

            if scheme_mlst.split('_', 1)[0] == scheme.split('_', 1)[0]:
                pass_qc = True
            else:
                if scheme == 'unknown' and scheme_mlst != '-':
                    pass_qc = True
                    logger.warning("Found {scheme_mlst} scheme for a species with unknown"
                                   " scheme".format(scheme_mlst=scheme_mlst))

                elif scheme == 'unknown' and scheme_mlst == '-':
                    pass_qc = True

                # a special warning was requested for yersinia
                elif species_genus == 'yersinia' and mlst_scheme_genus == 'yersinia':
                    pass_qc = True
                    logger.warning("Found a Yersinia scheme ({scheme_mlst}), but it is different from what it was"
                                   " expected ({scheme})".format(scheme_mlst=scheme_mlst, scheme=scheme))
                else:
                    if mlst_scheme_genus is not None and scheme_mlst == scheme == mlst_scheme_genus:
                        pass_qc = True
                    else:
                        logger.error("MLST scheme found ({scheme_mlst}) and provided ({scheme}) are not the"
                                     " same".format(scheme_mlst=scheme_mlst, scheme=scheme))
        else:
            pass_qc = True


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

        with open(".report.json", "w") as report:
            report.write(json.dumps(report_dic, separators=(",", ":")))

    else:
        logger.error("Sample {} did not run successfully".format(sample_id))

    # writing .status
    with open(".status", "w") as status:
        if pass_qc:
            status.write("pass")
        else:
            status.write("fail")


if __name__ == '__main__':

    main(SAMPLE_ID, ASSEMBLY, EXPECTED_SPECIES)
