diamond
=======

Purpose
-------

This component performs ``blastx`` or ``blastp`` with diamond. The database
used by diamond can be provided from the local disk or generated in the process.
This component uses the same output type as abricate with the same blast output
fields.

.. note::
    Software page: https://github.com/bbuchfink/diamond


Input/Output type
-----------------

- Input type: ``Fasta``
- Output type: None

.. note::
    The default input parameter for fasta data is ``--fasta``.

Parameters
----------

- ``pathToDb``: Provide full path for the diamond database. If none is provided
  then will try to fetch from the previous process. Default: None

- ``fastaToDb``: Provide the full path for the fasta to construct a diamond
  database. Default: None

- ``blastType``: Defines the type of blast that diamond will do. Can wither be
  blastx or blastp. Default: blastx

Published results
-----------------

- ``results/annotation/diamond*``: Stores the results of the abricate screening
  for each sample and for each specified database.

Published reports
-----------------

None.

Default directives
------------------

- ``diamond``:
    - ``container``: flowcraft/diamond
    - ``version``: 0.9.22-1
    - ``memory``: { 4.GB * task.attempt }
    - ``cpus``: 2