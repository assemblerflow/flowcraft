Integrity_coverage
==================

Purpose
-------

This module receives paired FastQ files, a genome size estimate and a
minimum coverage threshold and has three purposes while iterating over the
FastQ files:

- Checks the integrity of FastQ files (corrupted files).
- Guesses the encoding of FastQ files (this can be turned off in the opts argument).
- Estimates the coverage for each sample.


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

