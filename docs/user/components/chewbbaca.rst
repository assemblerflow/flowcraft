Chewbbaca
=========

Purpose
-------

This components runs the allele calling operation of ChewBBACA on a set
of fasta samples to perform a cg/wgMLST analysis

.. note::
    Software page: https://github.com/B-UMMI/chewBBACA

Input/Output type
------------------

- Input type: ``Fasta``
- Output type: None

.. note::
    The default input parameter for fasta data is ``--fasta``.

Parameters
----------

- ``chewbbacaQueue``: Specifiy a queue/partition for chewbbaca. This option
  is only used for grid schedulers.
- ``chewbbacaTraining``: Specify the full path to the prodigal training file
  of the corresponding species.
- ``schemaPath``: The path to the chewbbaca schema directory.
- ``schemaSelectedLoci``: The path to the selection of loci in the schema
  directory to be used. If not specified, all loci in the schema will be used.
- ``chewbbacaJson``: If set to True, chewbbaca's allele call output will be
  set to JSON format.
- ``chewbbacaToPhyloviz``: If set to True, the ExtractCgMLST module of
  chewbbaca will be executed after the allele calling.
- ``chewbbacaProfilePercentage``: Specifies the proportion of samples that
  must be present in a locus to save the profile.
- ``chewbbacaBatch``: Specifies whther a chewbbaca run will be performed on
  the complete input batch (all at the same time) or one by one.

Published results
-----------------

- ``results/chewbbaca_alleleCall``: The results of the allelecall for each
 sample.
- ``results/chewbbaca``: The cg/wgMLST schema prepared for phyloviz.

Published reports
-----------------

 None. 

Default directives
------------------

- ``chewbbaca``:
    - ``cpus``: 4
    - ``container``: mickaelsilva/chewbbaca_py3
    - ``version``: latest
- ``chewbbaca_batch``:
    - ``cpus``: 4
    - ``container``: mickaelsilva/chewbbaca_py3
    - ``version``: latest
- ``chewbbacaExtractMLST``:
    - ``container``: mickaelsilva/chewbbaca_py3
    - ``version``: latest
