seqsero2_reads
==============

Purpose
-------

Salmonella serotype determination from sequencing reads.

.. note::
    Software page: https://github.com/denglab/SeqSero2

Input/Output type
-----------------

- Input type: ``FastQ``
- Output type: None

.. note::
    The default input parameter for fastq data is ``--fastq``.

Published results
-----------------

- ``results/typing/seqsero2_reads``: Stores SeqSero2 results for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 1
- ``memory``: 1GB (dynamically increased with CPUs and on retry)
- ``container``: ummidock/seqsero2
- ``version``: alpha-test-1
- ``cache``: false
- ``scratch``: true

Advanced
--------

Reports JSON
^^^^^^^^^^^^

``typing``:
    - ``seqsero2_reads``: <serotype result>
