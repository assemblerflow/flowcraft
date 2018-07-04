remove_host
===========

Purpose
-------

This component performs a mapping procedure of FastQ files using a host
genome as referece (default: hg19). The procedure is carried out with
bowtie2 and samtools and aims to filter the reads that map to host genome.

.. note::
    - bowtie2 documentation can be found `here <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_.
    - samtools documentation can be found `here <http://www.htslib.org/doc/samtools-1.2.html>`_.

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: ``FastQ``

.. note::
    The default input parameter for fastq data is ``--fastq``.

Parameters
----------

- ``refIndex``: Specifies the reference indexes to be provided to bowtie2.
Default: '/index_hg19/hg19' (from docker image).


Published results
-----------------

- ``results/mapping/``: A `txt` file from bowtie2 with the mapping statistics.

Published reports
-----------------

None.

Default directives
------------------

- ``remove_host``:
    - ``cpus``: 3
    - ``memory``: 5GB (dynamically increased on retry)
    - ``container``: flowcraft/remove_host
    - ``version``: 2-0.1


Advanced
--------

Template
^^^^^^^^

:mod:`assemblerflow.templates.remove_host`
