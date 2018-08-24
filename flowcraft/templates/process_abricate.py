#!/usr/bin/env python3

"""
Purpose
-------

This module is intended parse the results of the Abricate for one or more
samples.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``abricate_files`` : Path to abricate output file.
    - e.g.: ``'abr_resfinder.tsv'``

Generated output
----------------

None


Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "26032018"
__template__ = "process_abricate-nf"

import re
import os
import json
import operator
import subprocess

from subprocess import PIPE

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


def __get_version_abricate():

    try:

        # Get abricate version
        cli = ["abricate", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        version = stdout.strip().split()[-1].decode("utf8")

    except Exception as e:
        logger.debug(e)
        version = "undefined"

    try:

        # Get abricate database versions
        cli = ["abricate", "--list"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        dbout, _ = p.communicate()

        databases = [[u.decode("utf8") for u in i.strip().split()]
                     for i in dbout.splitlines()][1:]

    except Exception as e:
        logger.debug(e)
        databases = "undefined"

    return {
        "program": "abricate",
        "version": version,
        "databases": databases
    }


if __file__.endswith(".command.sh"):
    ABRICATE_FILES = '$abricate_file'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("ABRICATE_FILE: {}".format(ABRICATE_FILES))


class Abricate:
    """Main parser for Abricate output files.

    This class parses one or more output files from Abricate, usually from
    different databases. In addition to the parsing methods, it also provides
    a flexible method to filter and re-format the content of the abricate
    files.

    Parameters
    ----------
    fls : list
       List of paths to Abricate output files.
    """

    def __init__(self, fls):

        self.storage = {}
        """
        dic: Main storage of Abricate's file content. Each entry corresponds
        to a single line and contains the keys::

            - ``log_file``: Name of the summary log file containing abricate
              results
            - ``infile``: Input file of Abricate.
            - ``reference``: Reference of the query sequence.
            - ``seq_range``: Range of the query sequence in the database
             sequence.
            - ``gene``: AMR gene name.
            - ``accession``: The genomic source of the sequence.
            - ``database``: The database the sequence came from.
            - ``coverage``: Proportion of gene covered.
            - ``identity``: Proportion of exact nucleotide matches.
        """

        self._key = 0
        """
        int: Arbitrary key for unique entries in the storage attribute
        """

        self.parse_files(fls)

    def parse_files(self, fls):
        """Public method for parsing abricate output files.

        This method is called at at class instantiation for the provided
        output files. Additional abricate output files can be added using
        this method after the class instantiation.

        Parameters
        ----------
        fls : list
            List of paths to Abricate files

        """

        for f in fls:
            # Make sure paths exists
            if os.path.exists(f):
                self._parser(f)
            else:
                logger.warning("File {} does not exist".format(f))

    def _parser(self, fl):
        """Parser for a single abricate output file.

        This parser will scan a single Abricate output file and populate
        the :py:attr:`Abricate.storage` attribute.

        Parameters
        ----------
        fl : str
            Path to abricate output file

        Notes
        -----
        This method will populate the :py:attr:`Abricate.storage` attribute
        with all compliant lines in the abricate output file. Entries are
        inserted using an arbitrary key that is set by the
        :py:attr:`Abricate._key` attribute.

        """

        with open(fl) as fh:

            for line in fh:
                # Skip header and comment lines
                if line.startswith("#") or line.strip() == "":
                    continue

                fields = line.strip().split("\t")

                try:
                    coverage = float(fields[8])
                except ValueError:
                    coverage = None
                try:
                    identity = float(fields[9])
                except ValueError:
                    identity = None

                try:
                    accession = fields[11]
                except IndexError:
                    accession = None

                self.storage[self._key] = {
                    "log_file": os.path.basename(fl),
                    "infile": fields[0],
                    "reference": fields[1],
                    "seq_range": (int(fields[2]), int(fields[3])),
                    "gene": fields[4],
                    "accession": accession,
                    "database": fields[10],
                    "coverage": coverage,
                    "identity": identity
                }

                self._key += 1

    @staticmethod
    def _test_truth(x, op, y):
        """ Test the truth of a comparison between x and y using an operator.

        If you want to compare '100 > 200', this method can be called as
        self._test_truth(100, ">", 200).

        Parameters
        ----------
        x : int
            Arbitrary value to compare in the left.
        op : str
            Comparison operator.
        y : int
            Arbitrary value to compare in the right.

        Returns
        -------
        x : bool
            The 'truthness' of the test.
        """

        ops = {
            ">": operator.gt,
            "<": operator.lt,
            ">=": operator.ge,
            "<=": operator.le,
            "==": operator.eq,
            "!=": operator.ne
        }

        return ops[op](x, y)

    def iter_filter(self, filters, databases=None, fields=None,
                    filter_behavior="and"):
        """General purpose filter iterator.

        This general filter iterator allows the filtering of entries based
        on one or more custom filters. These filters must contain
        an entry of the `storage` attribute, a comparison operator, and the
        test value. For example, to filter out entries with coverage below 80::

            my_filter = ["coverage", ">=", 80]

        Filters should always be provide as a list of lists::

            iter_filter([["coverage", ">=", 80]])
            # or
            my_filters = [["coverage", ">=", 80],
                          ["identity", ">=", 50]]

            iter_filter(my_filters)

        As a convenience, a list of the desired databases can be directly
        specified using the `database` argument, which will only report
        entries for the specified databases::

            iter_filter(my_filters, databases=["plasmidfinder"])

        By default, this method will yield the complete entry record. However,
        the returned filters can be specified using the `fields` option::

            iter_filter(my_filters, fields=["reference", "coverage"])

        Parameters
        ----------
        filters : list
            List of lists with the custom filter. Each list should have three
            elements. (1) the key from the entry to be compared; (2) the
            comparison operator; (3) the test value. Example:
                ``[["identity", ">", 80]]``.
        databases : list
            List of databases that should be reported.
        fields : list
            List of fields from each individual entry that are yielded.
        filter_behavior : str
            options: ``'and'`` ``'or'``
            Sets the behaviour of the filters, if multiple filters have been
            provided. By default it is set to ``'and'``, which means that an
            entry has to pass all filters. It can be set to ``'or'``, in which
            case one one of the filters has to pass.

        yields
        ------
        dic : dict
            Dictionary object containing a :py:attr:`Abricate.storage` entry
            that passed the filters.

        """

        if filter_behavior not in ["and", "or"]:
            raise ValueError("Filter behavior must be either 'and' or 'or'")

        for dic in self.storage.values():

            # This attribute will determine whether an entry will be yielded
            # or not
            _pass = False

            # Stores the flags with the test results for each filter
            # The results will be either True or False
            flag = []

            # Filter for databases
            if databases:
                # Skip entry if not in specified database
                if dic["database"] not in databases:
                    continue

            # Apply filters
            for f in filters:
                # Get value of current filter
                val = dic[f[0]]
                if not self._test_truth(val, f[1], f[2]):
                    flag.append(False)
                else:
                    flag.append(True)

            # Test whether the entry will pass based on the test results
            # and the filter behaviour
            if filter_behavior == "and":
                if all(flag):
                    _pass = True
            elif filter_behavior == "or":
                if any(flag):
                    _pass = True

            if _pass:
                if fields:
                    yield dict((x, y) for x, y in dic.items() if x in fields)
                else:
                    yield dic

    def get_filter(self, *args, **kwargs):
        """ Wrapper of the iter_filter method that returns a list with results

        It should be called exactly as in the `iter_filter`

        Returns
        -------
        _ : list
            List of dictionary entries that passed the filters in the
            `iter_filter` method.

        See Also
        --------
        iter_filter
        """

        return list(self.iter_filter(*args, **kwargs))


class AbricateReport(Abricate):
    """Report generator for single Abricate output files

    This class is intended to parse an Abricate output file from a single
    sample and database and generates a JSON report for the report webpage.

    Parameters
    ----------
    fls : list
       List of paths to Abricate output files.
    database : (optional) str
        Name of the database for the current report. If not provided, it will
        be inferred based on the first entry of the Abricate file.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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

    def get_plot_data(self):
        """ Generates the JSON report to plot the gene boxes

        Following the convention of the reports platform, this method returns
        a list of JSON/dict objects with the information about each entry in
        the abricate file. The information contained in this JSON is::

            {contig_id: <str>,
             seqRange: [<int>, <int>],
             gene: <str>,
             accession: <str>,
             coverage: <float>,
             identity: <float>
             }

        Note that the `seqRange` entry contains the position in the
        corresponding contig, not the absolute position in the whole assembly.

        Returns
        -------
        json_dic : list
            List of JSON/dict objects with the report data.
        """

        json_dic = {"plotData": []}
        sample_dic = {}
        sample_assembly_map = {}

        for entry in self.storage.values():

            sample_id = re.match("(.*)_abr", entry["log_file"]).groups()[0]
            if sample_id not in sample_dic:
                sample_dic[sample_id] = {}

            # Get contig ID using the same regex as in `assembly_report.py`
            # template
            contig_id = self._get_contig_id(entry["reference"])
            # Get database
            database = entry["database"]
            if database not in sample_dic[sample_id]:
                sample_dic[sample_id][database] = []

            # Update the sample-assembly correspondence dict
            if sample_id not in sample_assembly_map:
                sample_assembly_map[sample_id] = entry["infile"]

            sample_dic[sample_id][database].append(
                {"contig": contig_id,
                 "seqRange": entry["seq_range"],
                 "gene": entry["gene"].replace("'", ""),
                 "accession": entry["accession"],
                 "coverage": entry["coverage"],
                 "identity": entry["identity"],
                 },
            )

        for sample, data in sample_dic.items():
            json_dic["plotData"].append(
                {
                    "sample": sample,
                    "data": {"abricateXrange": data},
                    "assemblyFile": sample_assembly_map[sample]
                }
            )

        return json_dic

    def get_table_data(self):
        """

        Returns
        -------

        """

        gene_storage = {}
        json_dic = {"tableRow": []}
        logger.info("Generating JSON table data")

        # Collect the gene lists for each database
        for key, entry in self.storage.items():

            # Retrieve and initiate new sample entry, if not present already
            logger.debug("Retrieving sample if from: {}".format(
                entry["infile"]))
            sample_id = re.match("(.*)_abr", entry["log_file"]).groups()[0]
            database = entry["database"]

            if sample_id not in gene_storage:
                gene_storage[sample_id] = {}

            if database not in gene_storage[sample_id]:
                gene_storage[sample_id][database] = []

            gene_storage[sample_id][database].append(
                entry["gene"].replace("'", "").replace('"', '')
            )

        # For each database, create the JSON report
        for sample, table_data in gene_storage.items():

            json_dic["tableRow"].append({
                "sample": sample,
                "data": []
            })

            for db, gene_list in table_data.items():

                ind_json = {
                    "table": "abricate",
                    "header": db,
                    "value": len(gene_list),
                    "geneList": gene_list
                }
                json_dic["tableRow"][-1]["data"].append(ind_json)

        return json_dic

    def write_report_data(self):
        """Writes the JSON report to a json file
        """

        json_plot = self.get_plot_data()
        json_table = self.get_table_data()

        json_dic = {**json_plot, **json_table}

        with open(".report.json", "w") as json_report:
            json_report.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == '__main__':

    @MainWrapper
    def main(abr_file):

        abr = AbricateReport(fls=abr_file)
        abr.write_report_data()

    main(ABRICATE_FILES)
