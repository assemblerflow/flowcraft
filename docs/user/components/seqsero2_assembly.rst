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

- ``results/typing/seqsero2/assembly/seqsero2_assembly_ProcessID/``: Stores SeqSero2 results for each sample (for each seqsero2_assembly process).

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 1
- ``memory``: 1GB (dynamically increased on retry)
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
