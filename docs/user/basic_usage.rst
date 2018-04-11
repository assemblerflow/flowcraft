Basic Usage
===========

Assemblerflow has currently one execution mode, ``build``, that is used to
build the nextflow pipeline. However, more execution modes are being developed.

Building a pipeline
-------------------

Pipelines can be generated using the ``build`` execution mode of assemblerflow
and the ``-t`` option to specify the components inside quotes::

    assemblerflow build -t "trimmomatic fastqc spades" -o my_pipe.nf

This command will generate a linear pipeline with three components on the
current working directory. In addition to the main nextflow pipeline file
(``my_pipe.nf``), assemblerflow will write several auxiliary files that will
be necessary for the pipeline to run (templates, scripts, helpers, etc.).

