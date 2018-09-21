check_coverage
==============

Purpose
-------

This components estimates the coverage of a given sample based on the number
of base pairs in the FastQ files of a sample and on the expected genome size:

.. math::
    \frac{\text{number of base pairs}}{(\text{genome size} \times 1e^{6})}

If the estimated coverage of a given sample falls bellow the provided
minimum coverage threshold, the sample is filtered and does not proceed in the
pipeline.

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

Published results
-----------------

None.

Published reports
-----------------

- ``reports/coverage``: CSV table with estimated sequencing coverage for
  each sample.

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
    - ``Coverage``: Estimated coverage.
``fail``:
    - When estimated coverage is below the provided threshold.