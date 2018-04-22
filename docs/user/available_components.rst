.. _components:

Components
==========

These are the currently available assemblerflow components with a short
description of their tasks. For a more detailed information, follow the
links of each component.

Read Quality Control
--------------------

- :doc:`components/integrity_coverage`: Tests the integrity
  of the provided FastQ files, provides the option to filter FastQ files
  based on the expected assembly coverage and provides information about
  the maximum read length and sequence encoding.

- :doc:`components/fastqc`: Runs FastQC on paired-end FastQ files.

- :doc:`components/trimmomatic`: Runs Trimmomatic on paired-end FastQ files.

- :doc:`components/fastqc_trimmomatic`: Runs Trimmomatic on
  paired-end FastQ files informed by the FastQC report.

- :doc:`components/check_coverage`: Estimates the coverage for each sample and
  filters FastQ files according to a specified minimum coverage threshold.

Assembly
--------

- :doc:`components/spades`: Assembles paired-end FastQ files
  using SPAdes.

- :doc:`components/skesa`: Assembles paired-end FastQ files using
  skesa.

Post-assembly
-------------

- :doc:`components/process_spades`: Processes the assembly output
  from Spades and performs filtering base on quality criteria of GC content
  k-mer coverage and read length.

- :doc:`components/process_skesa`: Processes the assembly output
  from Skesa and performs filtering base on quality criteria of GC content
  k-mer coverage and read length.

- :doc:`components/assembly_mapping`: Performs a mapping
  procedure of FastQ files into a their assembly and performs filtering
  based on quality criteria of read coverage and genome size.

- :doc:`components/pilon`: Corrects and filters assemblies using Pilon.

Annotation
----------

- :doc:`components/prokka`: Performs assembly annotation using prokka.

- :doc:`components/abricate`: Performs anti-microbial gene screening using
  abricate.

MLST
----

- :doc:`components/mlst`: Checks the ST of an assembly using
  mlst.

- :doc:`components/chewbbaca`: Performs a cg/wgMLST analysis using ChewBBACA.


Reads typing
------------

- :doc:`components/seq_typing`: Determines the type of a given sample frm a set
  of reference sequences.
- :doc:`components/patho_typing`: *In silico* pathogenic typing from raw
  illumina reads.


Plasmids
--------

- :doc:`components/mapping_patlas`: Performs read mapping and generates a JSON
  input file for pATLAS.

- :doc:`components/mash_screen`: Performs mash screen against a reference index
  plasmid database and generates a JSON input file for pATLAS. This component
  searches for containment of a given sequence in read sequencing data.
  However if a different
  database is provided it can use mash screen for other purporses.

- :doc:`components/mash_dist`: Executes mash distance against a reference index
  plasmid database and generates a `JSON` for pATLAS. This component calculates
  pairwise distances between sequences (one from the database and the query
  sequence). However if a
  different database is provided it can use mash dist for other purposes.

