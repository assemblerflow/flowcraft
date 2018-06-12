trimmomatic
===========

Purpose
-------

This component runs Trimmomatic on paired-end FastQ files.

.. note::
    Software page: http://www.usadellab.org/cms/?page=trimmomatic

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

- ``adapters``: Provide a non-default fasta file containing the adapter
  sequences used to filter the FastQ files.
- ``trimSlidingWindow``: Perform sliding window trimming, cutting once the
  average quality within the window falls below a threshold.
- ``trimLeading``: Cut bases off the start of a read, if below a threshold
  quality.
- ``trimTrailing``: Cut bases of the end of a read, if below a threshold
  quality.
- ``trimMinLength``: Drop the read if it is below a specified length.

Published results
-----------------

- ``results/trimmomatic``: The trimmed FastQ files for each sample.

Published reports
-----------------

- ``reports/fastqc``: Stores the FastQC HTML reports for each sample.
- ``reports/fastqc/run_2/``: Stores the summary text files with the category
  results of FastQC for each sample.

Default directives
------------------

- ``cpus``: 2
- ``memory``: 4GB (dynamically increased on retry)
- ``container``: ummidock/trimmomatic
- ``version``: 0.36-2


Advanced
--------

Template
^^^^^^^^

:mod:`flowcraft.templates.trimmomatic`
:mod:`flowcraft.templates.trimmomatic_report`

Reports JSON
^^^^^^^^^^^^

``tableRow``:
    ``Trimmed (%)``: Percentage of trimmed nucleotides
``plotData``:
    ``sparkline``: Number of nucleotides after trimming
``badReads``: Number of discarded reads