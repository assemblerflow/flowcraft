FastQC
======

Purpose
-------

This components runs FastQC on paired-end FastQ files.

.. note::
    Software page: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

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
  sequences to screen overrepresented sequences against.

Published results
-----------------

None.

Published reports
-----------------

- ``reports/fastqc``: Stores the FastQC HTML reports for each sample.
- ``reports/fastqc/run_2/``: Stores the summary text files with the category
  results of FastQC for each sample.

Default directives
------------------

- ``cpus``: 2
- ``memory``: 4GB
- ``container``: ummidock/fastqc
- ``version``: 0.11.7-1

Advanced
--------

Template
^^^^^^^^

:mod:`assemblerflow.templates.fastqc_report`

Reports JSON
^^^^^^^^^^^^

``tableRow``:
    - ``Contigs``: Number of contigs
    - ``Assembled BP``: Number of assembled base pairs
``plotData``:
    - ``size_dist``: Distribution of contig size.
    - ``gcSliding``: Sliding window of the GC content along the genome
    - ``covSliding``: Sliding window of the coverage along the genome
