integrity_coverage
==================

Purpose
-------

This component is intended to test the integrity of the provided FastQ files.
It does so by attempting to parse uncompressed or compressed (``gz``, ``bz2``
or ``zip``) FastQ files (paired-end or single-end). During this parse, if the
FastQ files are not corrupt, it retrieves the following information:

- **sequence encoding**: Estimates the sequence encoding based on the quality
  scores. This information can then be passed to other components that might
  required it.
- **estimated coverage**: Provides a rough coverage estimation for each sample
  based on a user-provided genome size (see `Parameters`_). This estimation
  is essentially

  .. math::
      \frac{\text{number of base pairs}}{(\text{genome size} \times 1e^{6})}

  This information is written to the ``reports`` directory (See
  `Published reports`_)
- **maximum read length.**: Retrieves the maximum read length for each sample.

.. important::
    If the ``minCoverage`` parameter value is set to higher than 0, this
    component will filter samples with an estimated coverage below that
    threshold.

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: ``FastQ``

.. note::
    The default input parameter for FastQ data is ``--fastq``. You can change
    the ``--fastq`` parameter default pattern (``fastq/*_{1,2}.*``) according
    to input file names (e.g.: ``--fastq "path/to/fastq/*R{1,2}.*"``).

Parameters
----------

- ``genomeSize``: Genome size estimate for the samples. It is used to
  estimate the coverage and other assembly parameters and
  checks.
- ``minCoverage``: Minimum coverage for a sample to proceed. Can be set to
  0 to allow any coverage.

.. note::
    You can use these parameters as in the following example:
    ``--genomeSize 3``.

Published results
-----------------

None.

Published reports
-----------------

- ``reports/coverage``: CSV table with estimated sequencing coverage for
  each sample.
- ``reports/corrupted``: Text file with list of corrupted samples.

Default directives
------------------

None.


Advanced
--------

Template
^^^^^^^^

:mod:`flowcraft.templates.integrity_coverage`

Reports JSON
^^^^^^^^^^^^

``tableRow``:
    - ``Raw BP``: Number of nucleotides.
    - ``Reads``: Number of reads.
    - ``Coverage``: Estimated coverage.
``plotData``:
    - ``sparkline``: Number of nucleotides.
``warnings``:
    - When the enconding and/or phred score cannot be inferred from FastQ files.
``fail``:
    - When estimated coverage is below the provided threshold.