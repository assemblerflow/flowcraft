pilon
=====

Purpose
-------

This components Performs a mapping procedure of FastQ files into a their
assembly and performs filtering based on quality criteria of read coverage
and genome size.

.. note::
    Software page: https://github.com/broadinstitute/pilon

Input/Output type
------------------

- Input type: ``Fasta`` and ``FastQ``
- Output type: ``Fasta``

.. note::
    The default input parameter for fasta data is ``--fasta``.

Parameters
----------

None.

Published results
-----------------

- ``results/assembly/pilon``: Stores the polished fasta assemblies for each
  sample.

Published reports
-----------------

- ``reports/assembly/pilon``: Table with several summary statistics about the
  assembly for each sample.

Default directives
------------------

- ``pilon``:
    - ``cpus``: 4
    - ``memory``: 7GB (dynamically increased on retry)
    - ``container``: ummidock/pilon
    - ``version``: 1.22.0-2
- ``process_assembly_mapping``:
    - ``cpus``: 1
    - ``memory``: 7GB (dynamically increased on retry)
    - ``container``: ummidock/pilon
    - ``version``: 1.22.0-2

Advanced
--------

Template
^^^^^^^^

:mod:`flowcraft.templates.assembly_report`

Reports JSON
^^^^^^^^^^^^
``tableRow``:
    - ``Contigs``: Number of contigs.
    - ``Assembled BP``: Number of assembled base pairs.
``plotData``:
    - ``size_dist``: Distribution of contig size.
    - ``sparkline``: Number of assembled base pairs.
    - ``genomeSliding``:
        - ``gcData``: Genome sliding window of GC content.
        - ``covData``: Genome sliding window of read coverage depth.
        - ``window``: Size of sliding window
        - ``xbars``: Position of contigs along the genome sliding window.
        - ``assemblyFile``: Name of the input assembly file.
``warnings``:
    - When the number of contigs exceeds a given threshold.
``fail``:
    - When the genome size is below 80% or above 150% of the expected genome size.

