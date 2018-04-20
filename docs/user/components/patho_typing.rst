Patho_typing
==========

Purpose
-------

Patho_typing is a software for *in silico* pathogenic typing
directly from raw Illumina reads.

.. note::
    Software page: https://github.com/B-UMMI/patho_typing

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: None

Parameters
----------

- ``species``: Species name. Must be the complete species name with genus
  and species, e.g.: 'Yersinia enterocolitica'.

Published results
-----------------

- ``results/pathotyping/<sample id>``: Stores the results of patho_typing in
  text and tabular format.

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 4
- ``memory``: 4GB
- ``container``: ummidock/patho_typing
- ``version``: 0.3.0-1

Advanced
--------

Reports JSON
^^^^^^^^^^^^

``typing``:
    - ``pathotyping``: <typing result>