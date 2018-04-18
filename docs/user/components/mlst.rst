MLST
====

Purpose
-------

Checks the ST of an assembly using mlst.

Input/Output type
------------------

- Input type: ``Fasta``
- Output type: None

.. note::
    The default input parameter for fasta data is ``--fasta``.

Parameters
----------

- ``mlstSpecies``: Specifiy the expected species for MLST.

Published results
-----------------

- ``results/annotation/mlst``: Stores the results of the ST for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``container``: ummidock/mlst
