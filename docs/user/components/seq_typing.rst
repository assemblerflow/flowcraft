seq_typing
==========

Purpose
-------

Seq_typing is a software that determines the type of a given sample using a
read mapping approach against a set of reference sequences. Sample's reads
are mapped to the given reference sequences and, based on the length of the
sequence covered and it's depth of coverage, seq_typing decides which reference
sequence is more likely to be present and returns the type associated with
such sequences.

.. note::
    Software page: https://github.com/B-UMMI/seq_typing

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: None

Parameters
----------

- ``referenceFileO``: Fasta file containing reference sequences. If more
  than one file is passed via the 'referenceFileH parameter, a reference
  sequence for each file will be determined.
- ``referenceFileH``: Fasta file containing reference sequences. If more
  than one file is passed via the 'referenceFileO parameter, a reference
  sequence for each file will be determined.

Published results
-----------------

- ``results/seqtyping/<sample id>``: Stores the results of seq_typing in
  text and tabular format.

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 4
- ``memory``: 4GB
- ``container``: ummidock/seq_typing
- ``version``: 0.1.0-1

Advanced
--------

Reports JSON
^^^^^^^^^^^^

``typing``:
    - ``seqtyping``: <typing result>