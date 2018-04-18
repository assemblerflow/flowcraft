Pilon
=====

Purpose
-------

This components Performs a mapping procedure of FastQ files into a their
assembly and performs filtering based on quality criteria of read coverage
and genome size.

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