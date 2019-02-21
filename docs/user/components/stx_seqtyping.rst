stx_seqtyping
=============

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

- ``results/typing/stx_seqtyping/stx_seqtyping_ProcessID/``: Stores `ecoli_stx_subtyping.py` results for each sample (for each stx_seqtyping process).

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 2
- ``memory``: 1GB (2GB with 2 CPUs) (dynamically increased with CPUs and on retry)

 - Example: if using 4 CPUs, the allocated memory in the first attempt will be 4GB (`aloc_mem = mem * cpu * attempt`)

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
