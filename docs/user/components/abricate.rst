Abricate
========

Purpose
-------

This components performs anti-microbial gene screening using abricate. It
includes the default databases plus the ``virulencefinder`` database.

Input/Output type
------------------

- Input type: ``Fasta`` and ``FastQ``
- Output type: ``Fasta``

.. note::
    The default input parameter for fasta data is ``--fasta``.

Parameters
----------

- ``abricateDatabases``: Specify the databases for abricate.

Published results
-----------------

- ``results/annotation/abricate``: Stores the results of the abricate screening
  for each sample and for each specified database.

Published reports
-----------------

None.

Default directives
------------------

- ``abricate``:
    - ``container``: ummidock/abricate
    - ``version``: 0.8.0-1
- ``process_assembly_mapping``:
    - ``container``: ummidock/abricate
    - ``version``: 0.8.0-1

Advanced
--------

Template
^^^^^^^^

:mod:`assemblerflow.templates.process_abricate`


Reports JSON
^^^^^^^^^^^^

``tableRow``:
    - ``<database>``: List of gene names
``plotData``:
    - ``<database>``:
        - ``contig``: Contig ID
        - ``seqRange``: Genomic range of the contig
        - ``gene``: Gene name
        - ``accession``: Accession number
        - ``coverage``: Coverage of the match
        - ``identity``: Identity of the match