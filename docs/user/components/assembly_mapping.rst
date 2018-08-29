assembly_mapping
================

Purpose
-------

This component performs a mapping procedure of FastQ files using their assembly
as reference. The procedure is carried out with bowtie2 and samtools and aims
to filter the assembly based on quality criteria of read coverage
and expected genome size.

.. note::
    - bowtie2 documentation can be found `here <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_.
    - samtools documentation can be found `here <http://www.htslib.org/doc/samtools-1.2.html>`_.

Input/Output type
------------------

- Input type: ``Fasta`` and ``FastQ``
- Output type: ``Fasta``

.. note::
    The default input parameter for fasta data is ``--fasta``.

Parameters
----------

- ``minAssemblyCoverage``: In auto, the default minimum coverage for each
  assembled contig is 1/3 of the assembly mean coverage or 10x, if the mean
  coverage is below 10x.
- ``AMaxContigs``: A warning is issues if the number of contigs is over
  this threshold.
- ``genomeSize``: Genome size estimate for the samples. It is used to check
  the ratio of contig number per genome MB.

Published results
-----------------

None.

Published reports
-----------------

None.

Default directives
------------------

- ``assembly_mapping``:
    - ``cpus``: 4
    - ``memory``: 5GB (dynamically increased on retry)
    - ``container``: ummidock/bowtie2_samtools
    - ``version``: 1.0.0-2
- ``process_assembly_mapping``:
    - ``cpus``: 1
    - ``memory``: 5GB (dynamically increased on retry)
    - ``container``: ummidock/bowtie2_samtools
    - ``version``: 1.0.0-2

Advanced
--------

Template
^^^^^^^^

:mod:`flowcraft.templates.process_assembly_mapping`

Reports JSON
^^^^^^^^^^^^

``plotData``:
    - ``sparkline``: Total number of base pairs.
``warnings``:
    - When the number of contigs exceeds a provided threshold.
``fail``:
    - When the genome size is below 80% or above 150% of the expected genome size.