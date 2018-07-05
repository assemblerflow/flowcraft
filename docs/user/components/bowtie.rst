bowtie
======

Purpose
-------

This component performs a mapping procedure of FastQ files with a given reference.
The procedure is carried out with Bowtie2.
The reference can a set of Bowtie2 index files or a Fasta file. In the latter, the
necessary index will be created with Bowtie2-build and passed through to Bowtie2.

.. note::
    - Bowtie2 documentation can be found `here <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_.
    - Software page: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: ``Bam``

.. note::
    The default input parameter for Fastq data is ``--fastq``.

Parameters
----------

- ``reference``: Specifies the reference genome to be provided to to bowtie2-build.
- ``index``: Specifies the reference indexes to be provided to bowtie2.

.. note::
    An ``index`` OR a ``reference`` fasta file must be provided

Published results
-----------------

- ``results/mapping/bowtie``: Stores the results of the mapping for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``bowtie_build``:
    - ``cpus``: 1
    - ``memory``: 5GB (dynamically increased on retry)
    - ``container``: flowcraft/bowtie2_samtools
    - ``version``: 1.0.0-1
- ``bowtie``:
    - ``cpus``: 4
    - ``memory``: 5GB (dynamically increased on retry)
    - ``container``:flowcraft/bowtie2_samtools
    - ``version``: 1.0.0-1
