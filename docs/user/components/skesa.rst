skesa
=====

Purpose
-------

This components assembles paired-end FastQ files using the Skesa assembler.

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

None.

Published results
-----------------

- ``results/assembly/skesa``: Stores the fasta assemblies for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 4
- ``memory``: 5GB (dynamically increased on retry)
- ``container``: flowcraft/skesa
- ``version``: 2.3.0-1
- ``scratch``: true

Advanced
--------

Template
^^^^^^^^

:mod:`flowcraft.templates.skesa`