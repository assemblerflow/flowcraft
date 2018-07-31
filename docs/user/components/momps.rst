momps
========

Purpose
-------

This component performs Multi-Locus Sequence Typing (MLST) on Legionella pneumophila
from reads and assemblies.

.. note::
    Software page: https://github.com/bioinfo-core-BGU/mompS

Input/Output type
------------------

- Input type: ``fasta``
- Output type: None

.. note::
    The default input parameter for fasta data is ``--fasta``. This process
    also requires FastQ reads provided via the ``--fastq`` parameter.

Parameters
----------

None.

Published results
-----------------

- ``results/typing/momps``: Stores TSV files with the ST and allelic profiles
  for each strain.

Published reports
-----------------

None.

Default directives
------------------

- ``momps``:
    - ``container``: flowcraft/momps
    - ``version``: 0.1.0-4

Advanced
--------

Reports JSON
^^^^^^^^^^^^

``typing``:
    - ``momps``: <typing result>