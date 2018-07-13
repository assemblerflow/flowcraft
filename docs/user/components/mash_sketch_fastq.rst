mash_sketch_fastq
=================

Purpose
-------

This component performs mash sketch for fastq input files.

.. note::
    - MASH documentation can be found `here <https://mash.readthedocs.io/en/latest/>`_.


Input/Output type
------------------

- Input type: ``fastq``
- Output type: ``msh``


Parameters
----------

- ``kmerSize``: Parameter to set the kmer size for hashing. Default: 21.
  Default: false.

- ``sketchSize``: Parameter to set the number of hashes per sketch.
  Default: 1000.

- ``minKmer``: Minimum copies of each k-mer required to pass noise filter for
  reads. Default: 1.

- ``genomeSize``: Genome size (raw bases or with K/M/G/T). If specified, will
  be used for p-value calculation instead of an estimated size from k-mer
  content. Default: *false*, meaning that it won't be used. If you want to use
  it pass a number to this parameter.


Published results
-----------------

None.


Published reports
-----------------

None.


Default directives
------------------

- ``mashSketchFastq``:
    - ``container``: flowcraft/mash-patlas
    - ``version``: 1.4.1
