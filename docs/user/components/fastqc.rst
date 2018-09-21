fastqc
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

:mod:`flowcraft.templates.fastqc_report`

Reports JSON
^^^^^^^^^^^^

``plotData``:
    - ``base_sequence_quality``: Per base sequence quality data
        - (This structure is repeated for the other entries)
        - ``status``: Status of the category (PASS, WARN, etc)
        - ``data``: Plot data
    - ``sequence_quality``: Per sequence quality data
    - ``base_gc_content``: GC content distribution
    - ``base_n_content``: Per base N content
    - ``sequence_length_dist``: Distribution of sequence read length
    - ``per_base_sequence_content``: Per base sequence content
``warnings``:
    - List of failures or warnings for some non-sensitive FastQC categories
``fail``:
    - Failure message when sensitive FastQC categories fail or do not pass.
