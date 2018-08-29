downsample_fastq
================

Purpose
-------

downsample_fastq uses seqtk to subsample fastq read data to a target coverage depth
if the estimated coverage is higher than the provided target depth. When
no subsample is required, it outputs the original FastQ files.

.. note::
    Software page: https://github.com/lh3/seqtk

Input/Output type
------------------

- Input type: ``fastq``
- Output type: ``fastq``

Parameters
----------

- ``genomeSize``: Genome size estimate for the samples. It is used to
  estimate the coverage.
- ``depth``: The target depth to which the reads should be subsampled.

Published results
-----------------

- ``results/sample_fastq``: Stores the subsampled FastQ files

Published reports
-----------------

None.

Default directives
------------------

- ``cpus``: 1
- ``memory``: 4GB
- ``container``: flowcraft/seqtk
- ``version``: 1.3.0-3

Advanced
--------

Reports JSON
^^^^^^^^^^^^

``tableRow``:
    - ``Coverage``: Estimated coverage.