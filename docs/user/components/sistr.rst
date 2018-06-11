sistr
=====

Purpose
-------

Sistr (Salmonella In Silico Typing Resource) is a software for Serovar
predictions from whole-genome sequence assemblies by determination
of antigen gene and cgMLST gene alleles using BLAST. Mash MinHash can also be
used for serovar prediction.

.. note::
    Software page: https://github.com/peterk87/sistr_cmd

Input/Output type
------------------

- Input type: ``Fasta``
- Output type: None

.. note::
    The default input parameter for fasta data is ``--fasta``.

Parameters
----------

None

Published results
-----------------

- ``results/typing/sistr``: Stores the results of sistr in a tab file

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 4
- ``memory``: 4GB
- ``container``: ummidock/sistr_cmd
- ``version``: 1.0.2
