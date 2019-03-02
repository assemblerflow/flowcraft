.. _components:

Components
==========

These are the currently available FlowCraft components with a short
description of their tasks. For a more detailed information, follow the
links of each component.


Download
--------

- :doc:`components/reads_download`: Downloads reads from the SRA/ENA public
  databases from a list of accessions.

- :doc:`components/fasterq_dump`: Downloads reads from the SRA public databases
  from a list of accessions, using ``fasterq-dump``.

Reads Quality Control
--------------------

- :doc:`components/check_coverage`: Estimates the coverage for each sample and
  filters FastQ files according to a specified minimum coverage threshold.

- :doc:`components/fastqc`: Runs FastQC on paired-end FastQ files.

- :doc:`components/fastqc_trimmomatic`: Runs Trimmomatic on
  paired-end FastQ files informed by the FastQC report.

- :doc:`components/filter_poly`:  Runs PrinSeq on paired-end
  FastQ files to remove low complexity sequences.

- :doc:`components/integrity_coverage`: Tests the integrity
  of the provided FastQ files, provides the option to filter FastQ files
  based on the expected assembly coverage and provides information about
  the maximum read length and sequence encoding.

- :doc:`components/trimmomatic`: Runs Trimmomatic on paired-end FastQ files.

- :doc:`components/downsample_fastq`: Subsamples fastq files up to a target coverage
  depth.


Assembly
--------

- :doc:`components/megahit`: Assembles metagenomic paired-end FastQ files
  using megahit.

- :doc:`components/metaspades`: Assembles metagenomic paired-end FastQ files
  using metaSPAdes.

- :doc:`components/skesa`: Assembles paired-end FastQ files using
  skesa.

- :doc:`components/spades`: Assembles paired-end FastQ files
  using SPAdes.

Post-assembly
-------------

- :doc:`components/pilon`: Corrects and filters assemblies using Pilon.

- :doc:`components/process_skesa`: Processes the assembly output
  from Skesa and performs filtering base on quality criteria of GC content
  k-mer coverage and read length.

- :doc:`components/process_spades`: Processes the assembly output
  from Spades and performs filtering base on quality criteria of GC content
  k-mer coverage and read length.

Binning
-------

- :doc:`components/maxbin2`: An automatic tool for binning metagenomic sequences

Annotation
----------

- :doc:`components/abricate`: Performs anti-microbial gene screening using
  abricate.

- :doc:`components/card_rgi`: Performs anti-microbial resistance gene screening using
  CARD rgi (with contigs as input).

- :doc:`components/prokka`: Performs assembly annotation using prokka.

Distance Estimation
-------------------

- :doc:`components/mash_dist`: Executes mash distance against a reference index
  plasmid database and generates a `JSON` for pATLAS. This component calculates
  pairwise distances between sequences (one from the database and the query
  sequence). However if a different database is provided it can use mash dist
  for other purposes.

- :doc:`components/mash_screen`: Performs mash screen against a reference index
  plasmid database and generates a JSON input file for pATLAS. This component
  searches for containment of a given sequence in read sequencing data.
  However if a different database is provided it can use mash screen for other
  purposes.

- :doc:`components/fast_ani`: Performs pairwise comparisons between fastas,
given a multifasta as input for fastANI. It will split the multifasta into
single fastas that will then be provided as a matrix. The output will be the
all pairwise comparisons that pass the minimum of 50 aligned sequences with a
default length of 200 bp.

- :doc:`components/mash_sketch_fasta`: Performs mash sketch for fasta files.

- :doc:`components/mash_sketch_fastq`: Performes mash sketch for fastq files.

Mapping
-------

- :doc:`components/assembly_mapping`: Performs a mapping
  procedure of FastQ files into a their assembly and performs filtering
  based on quality criteria of read coverage and genome size.

- :doc:`components/bowtie`: Align short paired-end sequencing reads to long reference sequences

- :doc:`components/mapping_patlas`: Performs read mapping and generates a JSON
  input file for pATLAS.

- :doc:`components/remove_host`: Performs read mapping with bowtie2
  against the target host genome (default hg19) and removes the mapping reads

- :doc:`components/retrieve_mapped`: Retrieves the mapped reads of a previous
  bowtie2 mapping process.

Taxonomic Profiling
---------------------

- :doc:`components/kraken`: Performs taxonomic identification with kraken on FastQ files
  (minikrakenDB2017 as default database)

- :doc:`components/kraken2`: Performs taxonomic identification with kraken2 on FastQ files
  (minikraken2_v1_8GB as default database)

- :doc:`components/midas_species`: Performs taxonomic identification on FastQ files at the
  species level with midas (requires database)

Typing
------

- :doc:`components/chewbbaca`: Performs a core-genome/whole-genome Multilocus
  Sequence Typing analysis on an assembly using ChewBBACA.

- :doc:`components/metamlst`: Checks the Sequence Type of metagenomic reads using
  Multilocus Sequence Typing.

- :doc:`components/mlst`: Checks the Sequence Type of an assembly using
  Multilocus Sequence Typing.

- :doc:`components/patho_typing`: *In silico* pathogenic typing from raw
  illumina reads.

- :doc:`components/seq_typing`: Determines the type of a given sample from a set
  of reference sequences.

- :doc:`components/sistr`: Serovar predictions from whole-genome sequence assemblies
  by determination of antigen gene and cgMLST gene alleles.

- :doc:`components/momps`: Multi-locus sequence typing for Legionella pneumophila
  from assemblies and reads.
