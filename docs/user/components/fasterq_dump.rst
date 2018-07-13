fasterq_dump
============

Purpose
-------

This component downloads reads from the SRA public databases from a
list of accessions. This component uses ``fasterq-dump`` from
`NCBI sra-tools <https://github.com/ncbi/sra-tools>`_. ``fasterq-dump``
increases the download speed in comparison from ``fastq-dump`` by
**multi-threading** the extraction of FASTQ from SRA-accessions.
The reads for each accession are then emitted through
the main output of this component to any other component (or components) that
receive FastQ data.

Input/Output type
------------------

- Input type: ``accessions``
- Output type: ``fastq``

.. note::
    The default input parameter for Accessions data is ``--accessions``.

Parameters
----------

- ``option_file``: This options enables the *option-file* parameter of
``fasterq-dump``, allowing parameters to be passed.
- ``compress_fastq``: This options allows users to disable the compression of
the fastq files resulting from this component. The default (``true``) behavior
compresses the fastq files to *fastq.gz*.

Published results
-----------------

- ``reads/<accession>``: Stores the reads for each provided accession.

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 1
- ``memory``: 1GB
- ``container``: flowcraft/sra-tools
- ``version``: 2.9.1-1
