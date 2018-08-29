reads_download
==============

Purpose
-------

This component downloads reads from the SRA/ENA public databases from a
list of accessions. First, it tries to use `aspera connect`_ to download
reads, if a valid aspera key is provided. Otherwise it uses curl, which is
substantially slower. The reads for each accession are then emitted through
the main output of this component to any other component (or components) that
receive FastQ data.

.. _aspera connect: http://asperasoft.com/download_connect/

Input/Output type
------------------

- Input type: ``accessions``
- Output type: ``fastq``

.. note::
    The default input parameter for Accessions data is ``--accessions``.

Parameters
----------

- ``asperaKey``: Downloads fastq accessions using Aspera Connect
  by providing the private-key file 'asperaweb_id_dsa.openssh' normally found
  in ~/.aspera/connect/etc/asperaweb_id_dsa.openssh after the installation.

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
- ``container``: flowcraft/getseqena
- ``version``: 0.4.0-2
