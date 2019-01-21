stx_seqtyping_reads
===================

Purpose
-------

Gets Escherichia coli stx subtypes from raw reads via seq_typing.

.. note::
    Software page: https://github.com/B-UMMI/seq_typing#e-coli-stx-subtyping

Input/Output type
-----------------

- Input type: ``FastQ``
- Output type: None

.. note::
    The default input parameter for fastq data is ``--fastq``.

Published results
-----------------

- ``results/typing/stx_seqtyping_reads``: Stores `ecoli_stx_subtyping.py` results for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 4
- ``memory``: 1GB (dynamically increased with CPUs and on retry)
- ``container``: ummidock/seq_typing
- ``version``: 2.2-01
- ``cache``: false
- ``scratch``: true

Advanced
--------

Reports JSON
^^^^^^^^^^^^

``typing``:
    - ``stx_seqtyping_reads``: <stx type result>
