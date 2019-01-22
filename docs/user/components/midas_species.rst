midas_species
=============

Purpose
-------

This component performs MIDAS to assign taxonomic labels fro species to short DNA
sequences, usually obtained through metagenomic studies.

.. note::
    Software page: https://github.com/snayfach/MIDAS

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: None

.. note::
    The default input parameter for fastq data is ``--fastq``.

Parameters
----------

- ``midasDB``: Specifies MIDAS database. Default: /MidasDB/midas_db_v1.2

Published results
-----------------

- ``results/taxonomy/midas``: Stores the results of the screening
  for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``container``: flowcraft/midas
- ``version``: 1.3.2-0.1
- ``memory``: 2.Gb*task.attempt
- ``cpus``: 3