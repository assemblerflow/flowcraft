Available processes
===================

These are the currently supported assemblerflow processes with a short
description of their tasks. For a more detailed information, follow the
links of each process.

Annotation
----------

- :mod:`~assemblerflow.templates.prokka`: Performs assembly annotation using
  prokka.

- :mod:`~assemblerflow.templates.abricate`: Performs anti-microbial gene
  screening using abricate.

Assembly
--------

- :mod:`~assemblerflow.templates.spades`: Assembles paired-end FastQ files
  using SPAdes.

- :mod:`~assemblerflow.templates.skesa`: Assembles paired-end FastQ files using
  skesa.

Distance estimation
-------------------

FastQ QC
--------

- :mod:`~assemblerflow.templates.integrity_coverage`: Tests the integrity
  of the provided FastQ files, provides the option to filter FastQ files
  based on the expected assembly coverage and provides information about
  the maximum read length and sequence encoding.

- :mod:`~assemblerflow.templates.fastqc`: Runs FastQC on paired-end FastQ
  files. It has the option of filtering FastQ files based on quality control
  checks.

- :mod:`~assemblerflow.templates.trimmomatic`: Runs Trimmomatic on paired-end
  FastQ files.

- :mod:`~assemblerflow.templates.fastqc_trimmomatic`: Runs Trimmomatic on
  paired-end FastQ files informed by FastQC metrics.

- :mod:`~assemblerflow.templates.check_coverage`: Filters FastQ files based
  on the expected assembly coverage.

Mapping
-------

MLST
----

- :mod:`~assemblerflow.templates.mlst`: Checks the ST of an assembly using
  mlst.

- :mod:`~assemblerflow.templates.chewbbaca`: Performs a cg/wgMLST analysis
  using ChewBBACA.

Post-assembly
-------------

- :mod:`~assemblerflow.templates.process_spades`: Processes the assembly output
  from spades and performs filtering base on quality criteria of GC content
  k-mer coverage and read length.

- :mod:`~assemblerflow.templates.assembly_mapping`: Performs a mapping
  procedure of FastQ files into a their assembly and performs filtering
  based on quality criteria of read coverage and genome size.

- :mod:`~assemblerflow.templates.pilon`: Corrects and filters assemblies
  using Pilon.
