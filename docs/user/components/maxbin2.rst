maxbin2
=======

Purpose
-------

This component is an automated binning algorithm to recover genomes from multiple metagenomic datasets

.. note::
    Software page: https://sourceforge.net/projects/maxbin2/

Input/Output type
------------------

- Input type: ``Fasta``  and ``FastQ``
- Output type: ``Fasta``

.. note::
    The default input parameter for fasta is ``--fasta``. This process also requires FastQ files.
    If the FastQ files are input to any upstream process, those will be provided to maxbin2 automatically,
    if not, they can be provided with the parameter ``--fastq``.

Parameters
----------

- ``min_contig_lenght``: Minimum contig length. Default: 1000

- ``max_iteration``: Maximum Expectation-Maximization algorithm iteration number. Default: 50

- ``prob_threshold``: Probability threshold for EM final classification. Default: 0.9

Published results
-----------------

- ``results/maxbin2/``: Stores the results of the binning in a folder
  for each sample.

Published reports
-----------------

None.

Default directives
------------------

- ``container``: flowcraft/maxbin2
- ``version``: 2.2.4-1
- ``cpus``: 4
- ``memory``: 8.GB (dynamically increased on retry)


Template
^^^^^^^^

:mod:`assemblerflow.templates.maxbin2`