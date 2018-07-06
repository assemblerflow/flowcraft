retrieve_mapped
===============

Purpose
-------

This component retrieves the mapping reads of a previous bowtie mapping process.
The procedure is carried out with samtools and aims to retrieve the reads that map to target reference.

.. note::
    - samtools documentation can be found `here <http://www.htslib.org/doc/samtools-1.2.html>`_.

Input/Output type
------------------

- Input type: ``bam``
- Output type: ``FastQ``

.. note::
    This process has the ``bowtie2`` process as a dependency.

Parameters
----------

None

Published results
-----------------

- ``results/mapping/retrieve_mapped``: Contains the resulting ``FastQ`` files.

Published reports
-----------------

None.

Default directives
------------------

- ``remove_host``:
    - ``cpus``: 2
    - ``memory``: 5GB (dynamically increased on retry)
    - ``container``: flowcraft/bowtie2_samtools
    - ``version``: 1.0.0-1

