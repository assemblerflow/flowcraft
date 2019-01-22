kraken2
=======

Purpose
-------

This component performs Kraken2 to assign taxonomic labels to short DNA
sequences, usually obtained through metagenomic studies.

.. note::
    Software page: https://ccb.jhu.edu/software/kraken2/

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: txt

.. note::
    The default input parameter for fastq data is ``--fastq``.

Parameters
----------

- ``kraken2DB``: Specifies kraken2 database. Default: minikraken2_v1_8GB (in path inside the
default container)

Published results
-----------------

- ``results/taxonomy/kraken2``: Stores the results of the screening
  for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``container``: flowcraft/kraken2
- ``version``: 2.0.7-1
- ``cpus``: 3
- ``memory``: 5GB (dynamically increased on retry)
