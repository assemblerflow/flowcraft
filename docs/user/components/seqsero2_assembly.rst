seqsero2_assembly
=================

Purpose
-------

Salmonella serotype determination from genome assemblies.

.. note::
    Software page: https://github.com/denglab/SeqSero2

Input/Output type
-----------------

- Input type: ``Fasta``
- Output type: None

.. note::
    The default input parameter for fasta data is ``--fasta``.

Published results
-----------------

- ``results/typing/seqsero2_assembly``: Stores SeqSero2 results for each sample.

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
    - ``seqsero2_assembly``: <serotype result>
