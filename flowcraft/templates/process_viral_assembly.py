#!/usr/bin/env python3

import os
import json
import operator
from itertools import groupby

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper



"""
Purpose
-------

This module is intended to process the output of assembly process from a single
sample from the program Spades or Megahit for the report component.
The main input is an fasta file produced by the assembler.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id``: Sample Identification string.
    - e.g.: ``'SampleA'``
- ``assembly``: fasta file from the assembler.
    - e.g.: ``'spades.fasta'``
-  ``orfSize``: minimum contig size to be considered a complete ORF

Generated output
----------------
- ``.report.jason``: Data structure for the report

Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "11.09.2018"
__template__ = "viral_assembly-nf"

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLY = '$assembly'
    MINSIZE = '$min_size'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("MINSIZE: {}".format(MINSIZE))


class Assembly:
    """Class that parses and filters a Fasta assembly file

    This class parses an assembly fasta file, collects a number
    of summary statistics and metadata from the contigs, filters
    contigs based on user-defined metrics and writes filtered assemblies
    and reports.

    Parameters
    ----------
    assembly_file : str
        Path to assembly file.
    min_contig_len : int
        Minimum contig length when applying the initial assembly filter.
    min_kmer_cov : int
        Minimum k-mer coverage when applying the initial assembly.
        filter.
    sample_id : str
        Name of the sample for the current assembly.
    """

    def __init__(self, assembly_file, min_contig_len, min_kmer_cov,
                 sample_id, min_size):

        self.contigs = {}
        """
        dict: Dictionary storing data for each contig.
        """

        self.filtered_ids = []
        """
        list: List of filtered contig_ids.
        """

        self.min_gc = 0.05
        """
        float: Sets the minimum GC content on a contig.
        """

        self.sample = sample_id
        """
        str: The name of the sample for the assembly.
        """

        self.nORFs = 0
        """
        int: number of complete ORFs in the assembly.
        """

        self.report = {}
        """
        dict: Will contain the filtering results for each contig.
        """

        self.filters = [
            ["length", ">=", min_contig_len],
            ["kmer_cov", ">=", min_kmer_cov]
        ]
        """
        list: Setting initial filters to check when parsing the assembly file.
        This can be later changed using the 'filter_contigs' method.
        """

        # Parse assembly and populate self.contigs
        self._parse_assembly(assembly_file)

        #Gets the number of ORFs
        self.getORFs(assembly_file, min_size)

    def getORFs(self, assembly, min_size):

        f_open = open(assembly, "rU")

        entry = (x[1] for x in groupby(f_open, lambda line: line[0] == ">"))

        ORF = 0

        for header in entry:
            seq = "".join(s.strip() for s in entry.__next__())
            if len(seq) >= int(min_size):
                ORF += 1

        self.nORFs = ORF


    @staticmethod
    def _parse_coverage(header_str):
        """Attempts to retrieve the coverage value from the header string.

        It splits the header by "_" and then screens the list backwards in
        search of the first float value. This will be interpreted as the
        coverage value. If it cannot find a float value, it returns None.
        This search methodology is based on the strings of assemblers
        like spades and skesa that put the mean kmer coverage for each
        contig in its corresponding fasta header.

        Parameters
        ----------
        header_str : str
            String

        Returns
        -------
        float or None
            The coverage value for the contig. None if it cannot find the
            value in the provide string.
        """

        cov = None
        for i in header_str.split("_")[::-1]:
            try:
                cov = float(i)
                break
            except ValueError:
                continue

        return cov

    def _parse_assembly(self, assembly_file):
        """Parse an assembly fasta file.

        This is a Fasta parsing method that populates the
        :py:attr:`~Assembly.contigs` attribute with data for each contig in the
        assembly.

        The insertion of data on the self.contigs is done by the
        :py:meth:`Assembly._populate_contigs` method, which also calculates
        GC content and proportions.

        Parameters
        ----------
        assembly_file : str
            Path to the assembly fasta file.

        """

        # Temporary storage of sequence data
        seq_temp = []
        # Id counter for contig that will serve as key in self.contigs
        contig_id = 0
        # Initialize kmer coverage and header
        cov, header = None, None

        with open(assembly_file) as fh:

            logger.debug("Starting iteration of assembly file: {}".format(
                assembly_file))
            for line in fh:
                # Skip empty lines
                if not line.strip():
                    continue
                else:
                    # Remove whitespace surrounding line for further processing
                    line = line.strip()

                if line.startswith(">"):
                    # If a sequence has already been populated, save the
                    # previous contig information
                    if seq_temp:
                        # Use join() to convert string list into the full
                        # contig string. This is generally much more efficient
                        # than successively concatenating strings.
                        seq = "".join(seq_temp)

                        logger.debug("Populating contig with contig_id '{}', "
                                     "header '{}' and cov '{}'".format(
                                        contig_id, header, cov))
                        self._populate_contigs(contig_id, header, cov, seq)

                        # Reset temporary sequence storage
                        seq_temp = []
                        contig_id += 1

                    header = line[1:]
                    cov = self._parse_coverage(line)

                else:
                    seq_temp.append(line)

            # Populate last contig entry
            logger.debug("Populating contig with contig_id '{}', "
                         "header '{}' and cov '{}'".format(
                            contig_id, header, cov))
            seq = "".join(seq_temp)
            self._populate_contigs(contig_id, header, cov, seq)

    def _populate_contigs(self, contig_id, header, cov, sequence):
        """ Inserts data from a single contig into\
         :py:attr:`~Assembly.contigs`.

        By providing a contig id, the original header, the coverage that
        is parsed from the header and the sequence, this method will
        populate the :py:attr:`~Assembly.contigs` attribute.

        Parameters
        ----------
        contig_id : int
            Arbitrary unique contig identifier.
        header : str
            Original header of the current contig.
        cov : float
            The contig coverage, parsed from the fasta header
        sequence : str
            The complete sequence of the contig.

        """

        # Get AT/GC/N counts and proportions.
        # Note that self._get_gc_content returns a dictionary with the
        # information on the GC/AT/N counts and proportions. This makes it
        # much easier to add to the contigs attribute using the ** notation.
        gc_kwargs = self._get_gc_content(sequence, len(sequence))
        logger.debug("Populate GC content with: {}".format(gc_kwargs))

        self.contigs[contig_id] = {
            "header": header,
            "sequence": sequence,
            "length": len(sequence),
            "kmer_cov": cov,
            **gc_kwargs
        }

    @staticmethod
    def _get_gc_content(sequence, length):
        """Get GC content and proportions.

        Parameters
        ----------
        sequence : str
            The complete sequence of the contig.
        length : int
            The length of the sequence contig.

        Returns
        -------
        x : dict
            Dictionary with the at/gc/n counts and proportions

        """

        # Get AT/GC/N counts
        at = sum(map(sequence.count, ["A", "T"]))
        gc = sum(map(sequence.count, ["G", "C"]))
        n = length - (at + gc)

        # Get AT/GC/N proportions
        at_prop = at / length
        gc_prop = gc / length
        n_prop = n / length

        return {"at": at, "gc": gc, "n": n,
                "at_prop": at_prop, "gc_prop": gc_prop, "n_prop": n_prop}

    @staticmethod
    def _test_truth(x, op, y):
        """ Test the truth of a comparisong between x and y using an \
        ``operator``.

        If you want to compare '100 > 200', this method can be called as::

            self._test_truth(100, ">", 200).

        Parameters
        ----------
        x : int
            Arbitrary value to compare in the left
        op : str
            Comparison operator
        y : int
            Arbitrary value to compare in the rigth

        Returns
        -------
        x : bool
            The 'truthness' of the test
        """

        ops = {
            ">": operator.gt,
            "<": operator.lt,
            ">=": operator.ge,
            "<=": operator.le,
        }

        return ops[op](x, y)

    def filter_contigs(self, *comparisons):
        """Filters the contigs of the assembly according to user provided\
        comparisons.

        The comparisons must be a list of three elements with the
        :py:attr:`~Assembly.contigs` key, operator and test value. For
        example, to filter contigs with a minimum length of 250, a comparison
        would be::

            self.filter_contigs(["length", ">=", 250])

        The filtered contig ids will be stored in the
        :py:attr:`~Assembly.filtered_ids` list.

        The result of the test for all contigs will be stored in the
        :py:attr:`~Assembly.report` dictionary.

        Parameters
        ----------
        comparisons : list
            List with contig key, operator and value to test.

        """

        # Reset list of filtered ids
        self.filtered_ids = []
        self.report = {}

        gc_filters = [
            ["gc_prop", ">=", self.min_gc],
            ["gc_prop", "<=", 1 - self.min_gc]
        ]

        self.filters = list(comparisons) + gc_filters

        logger.debug("Filtering contigs using filters: {}".format(
            self.filters))

        for contig_id, contig in self.contigs.items():
            for key, op, value in list(comparisons) + gc_filters:
                if not self._test_truth(contig[key], op, value):
                    self.filtered_ids.append(contig_id)
                    self.report[contig_id] = "{}/{}/{}".format(key,
                                                               contig[key],
                                                               value)
                    break
                else:
                    self.report[contig_id] = "pass"

    def get_assembly_length(self):
        """Returns the length of the assembly, without the filtered contigs.

        Returns
        -------
        x : int
            Total length of the assembly.

        """

        return sum(
            [vals["length"] for contig_id, vals in self.contigs.items()
             if contig_id not in self.filtered_ids])

    def write_assembly(self, output_file, filtered=True):
        """Writes the assembly to a new file.

        The ``filtered`` option controls whether the new assembly will be
        filtered or not.

        Parameters
        ----------
        output_file : str
            Name of the output assembly file.
        filtered : bool
            If ``True``, does not include filtered ids.
        """

        logger.debug("Writing the filtered assembly into: {}".format(
            output_file))
        with open(output_file, "w") as fh:

            for contig_id, contig in self.contigs.items():
                if contig_id not in self.filtered_ids and filtered:
                    fh.write(">{}_{}\\n{}\\n".format(self.sample,
                                                     contig["header"],
                                                     contig["sequence"]))

    def write_report(self, output_file):
        """Writes a report with the test results for the current assembly

        Parameters
        ----------
        output_file : str
            Name of the output assembly file.

        """

        logger.debug("Writing the assembly report into: {}".format(
            output_file))
        with open(output_file, "w") as fh:

            for contig_id, vals in self.report.items():
                fh.write("{}, {}\\n".format(contig_id, vals))



@MainWrapper
def main(sample_id, assembly_file, minsize):
    """Main executor of the process_mapping template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    assembly: str
        Path to the fatsa file generated by the assembler.
    minsize: str
        Min contig size to be considered a complete ORF

    """

    logger.info("Starting assembly file processing")
    warnings = []
    fails = ""

    # Parse the spades assembly file and perform the first filtering.
    logger.info("Starting assembly parsing")
    assembly_obj = Assembly(assembly_file, 0, 0,
                            sample_id, minsize)

    if 'spades' in assembly_file:
        assembler = "SPAdes"
    else:
        assembler = "MEGAHIT"

    with open(".warnings", "w") as warn_fh:

        t_80 = int(minsize) * 0.8
        t_150 = int(minsize) * 1.5
        # Check if assembly size of the first assembly is lower than 80% of the
        # estimated genome size - DENV ORF has min 10k nt. If True, redo the filtering without the
        # k-mer coverage filter
        assembly_len = assembly_obj.get_assembly_length()
        logger.debug("Checking assembly length: {}".format(assembly_len))

        if assembly_len < t_80:

            logger.warning("Assembly size ({}) smaller than the minimum "
                           "threshold of 80% of expected genome size. "
                           "Applying contig filters without the k-mer "
                           "coverage filter".format(assembly_len))

            assembly_len = assembly_obj.get_assembly_length()
            logger.debug("Checking updated assembly length: "
                         "{}".format(assembly_len))
            if assembly_len < t_80:

                warn_msg = "Assembly size smaller than the minimum" \
                           " threshold of 80% of expected genome size: {}".format(
                                assembly_len)
                logger.warning(warn_msg)
                warn_fh.write(warn_msg)
                fails = warn_msg

        if assembly_len > t_150:

            warn_msg = "Assembly size ({}) larger than the maximum" \
                       " threshold of 150% of expected genome size.".format(
                            assembly_len)
            logger.warning(warn_msg)
            warn_fh.write(warn_msg)
            fails = warn_msg


    # Write json report
    with open(".report.json", "w") as json_report:
        json_dic = {
            "tableRow": [{
                "sample": sample_id,
                "data": [
                    {"header": "Contigs ({})".format(assembler),
                     "value": len(assembly_obj.contigs),
                     "table": "assembly",
                     "columnBar": True},
                    {"header": "Assembled BP ({})".format(assembler),
                     "value": assembly_len,
                     "table": "assembly",
                     "columnBar": True},
                    {"header": "ORFs",
                     "value": assembly_obj.nORFs,
                     "table": "assembly",
                     "columnBar":False}
                ]
            }],
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

    with open(".status", "w") as status_fh:
        status_fh.write("pass")



if __name__ == '__main__':

    main(SAMPLE_ID, ASSEMBLY, MINSIZE)

