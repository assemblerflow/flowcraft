metamlst
========

Purpose
-------

Checks the ST of metagenomic reads using mlst.

.. note::
    Software page: https://bitbucket.org/CibioCM/metamlst

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: None

.. note::
    The default input parameter for fastq data is ``--fastq``.

Parameters
----------

- ``metamlstDB``: Specifiy the metamlst database (full path) for MLST checking

- ``metamlstDB_index``: Specifiy the Bowtie2 metamlst database index (full path) for MLST checking

Published results
-----------------

- ``results/annotation/metamlst``: Stores the results of the ST for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``container``: flowcraft/metamlst
- ``version``: 1.1-1
- ``memory``: 4.Gb * task.attempt

