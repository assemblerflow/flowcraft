..    include:: <isonum.txt>

Overview
========

FlowCraft is an assembler of pipelines written in  nextflow_ for
analyses of genomic data. The premisse is simple:

Software are container blocks |rarr| Build your lego-like pipeline |rarr| Execute it (almost) anywhere.

What is Nextflow
::::::::::::::::

If you do not know nextflow, be sure to check it out. It's an awesome
framework based on the dataflow programming model used for building
parallelized, scalable and reproducible workflows using software containers.
It provides an abstraction layer between the execution and the logic of the
pipeline, which means that the same pipeline code can be executed on
multiple platforms, from a local laptop to clusters managed with SLURM, SGE,
etc. These are quite attractive features since genomic pipelines are
increasingly executed on large computer clusters to handle large volumes
of data and/or tasks. Moreover, portability and reproducibility are becoming
central pillars in modern data science.

What FlowCraft does
:::::::::::::::::::

FlowCraft is a python engine that automatically builds nextflow pipelines
by assembling pre-made ready-to-use :ref:`components <components>`. These components are modular
pieces of software or scripts, such as ``fastqc``, ``trimmomatic``, ``spades``,
etc, that are written for nextflow and have a set of attributes, such as
input and output types, parameters, directives, etc. This modular nature
allows them to be freely connected as long as they respect some basic rules,
such as the input type of a component must match with the output type of
the preceding component. In this way, nextflow processes can be
written only once, and FlowCraft is the magic glue that connects them,
handling the linking and forking of channels automatically. Moreover, each
component is associated with a docker image, which means that there is no
need to install any dependencies at all and all software runs on a
transparent and reliable box. To illustrate:

- A linear genome assembly pipeline can be easily built using FlowCraft
  with the following pipeline string::

    trimmomatic fastqc spades

Which will generate all the necessary files to run the nextflow
pipeline on any linux system that has nextflow and a container engine.

- You can easily add more components to perform assembly polishing, in this
  case, ``pilon``::

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
  accessions (downloaded from public databases) in a single workflow::

    download_reads trimmomatic={'extra_input':'reads'} fastqc skesa pilon

This pipeline can be executed by providing a file with accession numbers
(``--accessions`` parameter by default) **and** fastq reads, using the
``--reads`` parameter defined with the ``extra_input`` directive.


Who is FlowCraft for
::::::::::::::::::::

FlowCraft can be useful for bioinformaticians with varied levels of expertise
that need to executed genomic pipelines often and potentially in different
platforms. Building and executing pipelines requires no programming knowledge,
but familiarization with nextflow is highly recommended to take full advantage
of the generated pipelines.

At the moment, the available pre-made processes are mainly focused on
bacterial genome assembly simply because that was how we started.
However, our goal is to expand the library of existing components to other
commonly used tools in the field of genomics and to widen the applicability
and usefulness of FlowCraft pipelines.

Why not just write a Nextflow pipeline?
:::::::::::::::::::::::::::::::::::::::

In many cases, building a static nextflow pipeline is sufficient for our goals.
However, when building our own pipelines, we often felt the need to add
dynamism to this process, particularly if we take into account how fast new
tools arise and existing ones change. Our biological goals also change over
time and we might need different pipelines to answer different questions.
FlowCraft makes this very easy by having a set of pre-made and ready-to-use
components that can be freely assembled. By using components (``fastqc``,
``trimmomatic``) as its atomic elements, very complex pielines that take
full advantage of nextflow can be built with little effort. Moreover,
these components have explicit and standardized
input and output types, which means that the addition of new modules does not
require any changes in the existing code base. They just need to take into
account how data will be received by the process and how data may be emitted
from the process, to ensure that it can link with other components.

**However, why not both?**

FlowCraft generates a complete Nextflow pipeline file, which ca be used
as a starting point for your customized processes!

.. _nextflow: https://www.nextflow.io/