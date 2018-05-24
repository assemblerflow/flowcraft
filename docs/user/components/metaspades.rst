metaSPAdes
==========

Purpose
-------

This components assembles metagenomic paired-end FastQ files using the metaSPAdes assembler.

.. note::
    Software page: http://bioinf.spbau.ru/spades

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

- ``metaspadesKmers``: If 'auto' the metaSPAdes k-mer lengths will be determined
  from the maximum read length of each assembly. If 'default', metaSPAdes will
  use the default k-mer lengths.

Published results
-----------------

- ``results/assembly/metaspades``: Stores the fasta assemblies for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 4
- ``memory``: 5GB (dynamically increased on retry)
- ``container``: ummidock/spades
- ``version``: 3.11.1-1
- ``scratch``: true

Advanced
--------

Template
^^^^^^^^

:mod:`assemblerflow.templates.metaspades`