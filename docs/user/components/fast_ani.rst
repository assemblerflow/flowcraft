fast_ani
========

Purpose
-------

This component performs pairwise comparisons between fastas,
given a multifasta as input for fastANI. It will split the multifasta into
single fastas that will then be provided as a matrix. The output will be the
all pairwise comparisons that pass the minimum of 50 aligned sequences with a
default length of 200 bp.

Input/Output type
------------------

- Input type: ``fasta``
- Output type: ``None``


Parameters
----------

- ``fragLen``: Sets the minimum size of the fragment to be passed to
`--fragLen` argument of fastANI.


Published results
-----------------

- ``results/fast_ani/``: A text file with the extension `.out`, which has all
the pairwise comparisons between sequences, reporting ANI.


Published reports
-----------------

None.


Default directives
------------------

- ``fastAniMatrix``:
    - ``container``: flowcraft/fast_ani
    - ``version``: 1.1.0-2
    - ``cpus``: 20
    - ``memory``: { 30.GB * task.attempt }