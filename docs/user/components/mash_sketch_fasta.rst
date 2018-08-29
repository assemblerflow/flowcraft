mash_sketch_fasta
=================

Purpose
-------

This component performs mash sketch for fasta input files.

.. note::
    - MASH documentation can be found `here <https://mash.readthedocs.io/en/latest/>`_.


Input/Output type
------------------

- Input type: ``fasta``
- Output type: ``msh``


Parameters
----------

- ``kmerSize``: Parameter to set the kmer size for hashing. Default: 21.
  Default: false.

- ``sketchSize``: Parameter to set the number of hashes per sketch.
  Default: 1000.


Published results
-----------------

None.


Published reports
-----------------

None.


Default directives
------------------

- ``mashSketchFasta``:
    - ``container``: flowcraft/mash-patlas
    - ``version``: 1.4.1
