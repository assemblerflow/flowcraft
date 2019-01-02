kraken
======

Purpose
-------

This component performs Kraken to assign taxonomic labels to short DNA
sequences, usually obtained through metagenomic studies.

.. note::
    Software page: https://ccb.jhu.edu/software/kraken/

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: None

.. note::
    The default input parameter for fastq data is ``--fastq``.

Parameters
----------

- ``krakenDB``: Specifies kraken database. Default: minikraken_20171013_4GB (in path)

Published results
-----------------

- ``results/taxonomy/kraken``: Stores the results of the screening
  for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``container``: flowcraft/kraken
- ``version``: 1.0-0.1