Overview
========

Assemblerflow is an assembler of pipelines written in  nextflow_ for several
analyses of genomic data.

What is Nextflow
::::::::::::::::

If you do not know nextflow, be sure to check it out. It's an awesome
framework based on the dataflow programming model used for building
parallelized, scalable and reproducible workflows using software containers.
It provides an abstraction layer between the execution and the logic of the
pipeline, which means that the same pipeline code can be executed in
multiple platforms, such as SLURM, SGE, etc. These are quite attractive features
since genomic pipelines are increasingly executed on large computer clusters
to handle large volumes of data and/or tasks. Moreover, portability and
reproducibility are becoming central pillars in modern data science.

What Assemblerflow does
:::::::::::::::::::::::

Assemblerflow is a python engine that automatically builds nextflow pipelines
by assembling pre-made ready-to-use components. These components are modular
pieces of software or scripts, such as ``fastqc``, ``trimmomatic``, ``spades``,
etc, that are written for nextflow and have a set of attributes, such as
input and output types, parameters, directives, etc. This modular nature
allows them to be freely connected as long as they respect some basic rules,
such as the input type of a component must match with the output type of
the preceding component. In this way, portions of nextflow pipelines can be
written only once, and then freely used to create any number of pipeline
variations. Moreover, each component is associated with a docker image,
which means that there is no need to install any dependencies at all and
all software runs on a transparent and reliable box.
To illustrate:

- A linear genome assembly pipeline can be easily built using assemblerflow
  with the following pipeline string::

    trimmomatic fastqc spades

Which will generate all the necessary files to run the nextflow
pipeline on any linux system that has nextflow and a container engine.

- You can easily add more components to perform assembly polishing::

    trimmomatic fastqc spades pilon

- If a new assembler comes along and you want to switch that component in the
  pipeline, its as easy as replacing ``spades`` (or any other component)::

    trimmomatic fastqc skesa pilon

- And you can also fork the output of a component into multiple ones. For
  instance, we could annotate the resulting assemblies with multiple software::

    trimmomatic fastqc spades pilon (abricate | prokka)

- Or fork the execution of a pipeline early on to compare different software::

    trimmomatic fastqc (spades pilon | skesa pilon)

This will fork the output of ``fastqc`` into ``spades`` and ``skesa``, and
the pipeline will proceed independently in these two new 'lanes'.

- Directives for each process can be dynamically set when building the pipeline,
  such as the cpu/RAM usage or the software version::

    trimmomatic={'cpus':'4'} fastqc={'version':'0.11.5'} skesa={'memory':'10GB'} pilon (abricate | prokka)

- And extra input can be directly inserted in any part of the pipeline. For
  example, it is possible to assemble genomes from both fastq files and SRR
  accessions (downloaded from public databases)::

    download_reads trimmomatic={'extra_input':'reads'} fastqc skesa pilon

This pipeline can be executed by providing a file with accession numbers
(``--accessions`` parameter by default) **and** fastq reads, using the
``--reads`` parameter.


Who is Assemblerflow for
::::::::::::::::::::::::

Assemblerflow can be useful for bioinformaticians of varied levels of expertise
that need to executed genomic pipelines often and potentially in different
platforms. Building and executing pipelines requires no programming knowledge
but familiarization with nextflow is highly recommended to take full advantage
of the generated pipelines.

At the moment, the available pre-made processes are mainly focused on
bacterial genome assembly simply because that was how we started.
However, our goal is to expand the library of existing components to other
commonly used tools in the field of genomics and to widen the applicability
and usefulness of assemblerflow pipelines.

.. _nextflow: https://www.nextflow.io/