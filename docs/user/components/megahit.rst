megahit
=======

Purpose
-------

This components assembles metagenomic paired-end FastQ files using the megahit assembler.

.. note::
    Software page: https://github.com/voutcn/megahit

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: ``Fasta``

.. note::
    The default input parameter for FastQ data is ``--fastq``. You can change
    the ``--fastq`` parameter default pattern (``fastq/*_{1,2}.*``) according
    to input file names (e.g.: ``--fastq "path/to/fastq/*R{1,2}.*"``).

Parameters
----------

- ``megahitKmers``: If 'auto' the megahit k-mer lengths will be determined
  from the maximum read length of each assembly. If 'default', megahit will
  use the default k-mer lengths.

- ``fastg``: When true, it converts megahit intermediate contigs into fastg.
  Default: False



Published results
-----------------

- ``results/assembly/megahit``: Stores the fasta assemblies for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 4
- ``memory``: 5GB (dynamically increased on retry)
- ``container``: cimendes/megahit
- ``version``: v1.1.3-0.1
- ``scratch``: true

Advanced
--------

Template
^^^^^^^^

:mod:`assemblerflow.templates.megahit`