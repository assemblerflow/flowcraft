filter_poly
===========

Purpose
-------

This component removes low complexity sequence from read data
using PrinSeq.

.. note::
    Software page: http://prinseq.sourceforge.net/

Input/Output type
------------------

- Input type: ``FastQ``
- Output type: ``FastQ`

.. note::
    The default input parameter for fastq data is ``--fastq``.

Parameters
----------

- ``adapter``: Pattern to filter the reads. Please separate parameter values with a space
    and separate new parameter sets with semicolon (;). Parameters are defined by two values:
    the pattern (any combination of the letters ATCGN), and the number of repeats or percentage
    of occurence. Default: A 50%; T 50%; N 50%

Published results
-----------------

None.

Published reports
-----------------

None.

Default directives
------------------

- ``container``: flowcraft/prinseq
- ``version``: 0.20.4-1
- ``memory``: 4.GB * task.attempt
- ``cpus``: 1


