process_skesa
==============

Purpose
-------

This components processes the assembly resulting from the Skesa software and,
optionally, filters contigs based on user-provide parameters.

Input/Output type
------------------

- Input type: ``Fasta``
- Output type: ``Fasta``

.. note::
    The default input parameter for fasta data is ``--fasta``.

Parameters
----------

- ``skesaMinKmerCoverage``: Minimum contigs K-mer coverage. After assembly
  only keep contigs with reported k-mer coverage equal or above this value.
- ``skesaMinContigLen``: Filter contigs for length greater or equal than
  this value.
- ``skesaMaxContigs``: Maximum number of contigs per 1.5 Mb of expected
  genome size.

Published results
-----------------

None.

Published reports
-----------------

- ``reports/assembly/skesa_filter``: The filter status for each contig and
  each sample. If any contig does not pass the filters, it reports which 
  filter type it failed and the corresponding value.

Default directives
------------------

- ``container``: ummidock/skesa
- ``version``: 0.2.0-3

Advanced
--------

Template
^^^^^^^^

:mod:`flowcraft.templates.process_assembly`

Reports JSON
^^^^^^^^^^^^

``tableRow``:
    - ``Contigs (<assembler>)``: Number of contigs.
    - ``Assembled BP (<assembler>)``: Number of assembled base pairs.
``warnings``:
    - When the number of contigs exceeds a given threshold.
``fail``:
    - When the genome size is below 80% or above 150% of the expected genome size.
