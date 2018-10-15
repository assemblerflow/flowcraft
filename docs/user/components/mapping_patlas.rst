mapping_patlas
==============

Purpose
-------

This component performs mapping (using `bowtie2` and `samtools`) against a
plasmid database in order to find
plasmids contained in high throughoput sequencing data. Then, the resulting file
can be imported into `pATLAS <http://www.patlas.site/>`_.

.. note::
    - pATLAs documentation can be found `here <https://tiagofilipe12.gitbooks.io/patlas/content/>`_.
    - bowtie2 documentation can be found `here <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_.
    - samtools documentation can be found `here <http://www.htslib.org/doc/samtools-1.2.html>`_.

Input/Output type
------------------

- Input type: ``fastq``
- Output type: ``json``


Parameters
----------

- ``max_k``: Sets the k parameter for bowtie2 allowing to make multiple mappings
  of the same read against several hits on the query sequence or sequences.
  Default: 10949.

- ``trim5``: Sets trim5 option for bowtie. This will become legacy with QC
  integration, but it enables to trim 5' end of reads to be mapped with bowtie2.
  Default: 0

- ``lengthJson``: A dictionary of all the lengths of reference sequences.
  Default: 'jsons/*_length.json' (from docker image).

- ``refIndex``: Specifies the reference indexes to be provided to bowtie2.
  Default: '/ngstools/data/indexes/bowtie2idx/bowtie2.idx' (from docker image).

- ``samtoolsIndex``: Specifies the reference indexes to be provided to samtools.
  Default: '/ngstools/data/indexes/fasta/samtools.fasta.fai' (from docker image).


Published results
-----------------

- ``results/mapping/``: A `JSON` file that can be imported to `pATLAS <http://www.patlas.site/>`_
  with the results from mapping.


Published reports
-----------------

None.


Default directives
------------------

- ``mappingBowtie``:
    - ``container``: flowcraft/mapping-patlas
    - ``version``: 1.6.0-1
- ``samtoolsView``:
    - ``container``: flowcraft/mapping-patlas
    - ``version``: 1.6.0-1
- ``jsonDumpingMapping``:
    - ``container``: flowcraft/mapping-patlas
    - ``version``: 1.6.0-1
