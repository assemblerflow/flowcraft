spades
======

Purpose
-------

This components assembles paired-end FastQ files using the Spades assembler.

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

- ``spadesMinCoverage``: The minimum number of reads to consider an edge in
  the de Bruijn graph during the assembly
- ``spadesMinKmerCoverage``: Minimum contigs K-mer coverage. After assembly
  only keep contigs with reported k-mer coverage equal or above this value
- ``spadesKmers``: If 'auto' the SPAdes k-mer lengths will be determined
  from the maximum read length of each assembly. If 'default', SPAdes will
  use the default k-mer lengths.

Published results
-----------------

- ``results/assembly/spades``: Stores the fasta assemblies for each sample.

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

:mod:`flowcraft.templates.spades`