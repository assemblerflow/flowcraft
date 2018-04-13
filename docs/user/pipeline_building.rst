Pipeline building
=================

Assemblerflow offers a few extra features when building pipelines using the
``build`` execution mode.

Raw input types
---------------

Forks
-----

The output of any component in an assemblerflow pipeline can be forked into
two or more components, using the following fork syntax::

    trimmomatic fastqc (spades | skesa)

.. image:: ../resources/fork_1.png
   :scale: 80 %
   :align: center

In this example, the output of ``fastqc`` will be fork into two new *lanes*,
which will proceed independently from each other. In this syntax, a fork is
triggered by the ``(`` symbol (and the corresponding closing ``)``) and each
lane will be separated by a ``|`` symbol. There is no limitation to the number
of forks or lanes that a pipeline has. For instance, we could add more
components after the ``skesa`` module, including another fork::

    trimmomatic fastqc (spades | skesa pilon (abricate | prokka | chewbbaca))

.. image:: ../resources/fork_2.png
   :scale: 80 %
   :align: center

In this example, data will be forked after ``fastqc`` into two new lanes,
processed by ``spades`` and ``skesa``. In the skesa lane, data will continue
to flow into the ``pilon`` component and its output will fork into three new
lanes.

It is also possible to start a fork at the beggining of the pipeline, which
basically means that the pipeline will have multiple starting points. If we
want to provide the raw input two multiple process, the fork syntax can start
at the beginning of the pipeline::

    (seq_typing | trimmomatic fastqc (spades | skesa))

.. image:: ../resources/fork_3.png
   :scale: 80 %
   :align: center

In this case, since both initial components (``seq_typing`` and
``integrity_coverage``) received fastq files as input, the data provided
via the ``--fastq`` parameter will be forked and provided to both processes.

.. note::
    Some components have dependencies which need to be included previously
    in the pipeline. For instance, ``trimmomatic`` requires
    ``integrity_coverage`` and ``pilon`` requires ``assembly_mapping``. By
    default, assemblerflow will insert any missing dependencies right before
    the process, which is why these components appear in the figures above.

.. warning::
    Pay special attention to the syntax of the pipeline string when using
    forks. However, when unable to parse it, assemblerflow will do its best
    to inform you where the parsing error occurred.

Directives
----------

Several directives with information on cpu usage, RAM, version, etc. can be
specified for each individual component when building the pipeline using the
``={}`` notation. These
directives are written to the ``resources.config`` and
``containers.config`` files that are generated in the pipeline directory. You
can pass any of the directives already supported by nextflow (https://www.nextflow.io/docs/latest/process.html#directives),
but the most commonly used include:

    - ``cpus``
    - ``memory``
    - ``queue``

In addition, you can also pass the ``container`` and ``version`` directives
which are parsed by assemblerflow to dynamically change the container and/or
version tag of any process.

Here is an example where we specify cpu usage, allocated memory and container
version in the pipeline string::

    assemblerflow build -t "fastqc={'version':'0.11.5'} \
                            trimmomatic={'cpus':'2'} \
                            spades={'memory':'\'10GB\''}" -o my_pipeline.nf

When a directive is not specified, it will assume the default value of the
nextflow directive.

.. warning::
    Take special care not to include any white space characters inside the
    directives field. Common mistakes occur when specifying directives like
    ``fastqc={'version': '0.11.5'}``.

.. note::
    The values specified in these directives are placed in the
    respective config files exactly as they are. For instance,
    ``spades={'memory':'10GB'}"`` will appear in the config as
    ``spades.memory = 10Gb``, which will raise an error in nextflow because
    ``10Gb`` should be a string. Therefore, if you want a string you'll need to add
    the ``'`` as in this example: ``spades={'memory':'\'10GB\''}"``. The
    reason why these directives are not automatically converted is to allow
    the specification of dynamic computing resources, such as
    ``spades={'memory':'{10.Gb*task.attempt}'}"``

Extra inputs
------------

By default, only the first process (or processes) in a pipeline will receive
the raw input data provided by the user. However, the ``extra_input`` special
directive allows one or more processes to receive input from an additional parameter
that is provided by the user::

    reads_download integrity_coverage={'extra_input':'local'} trimmomatic spades

The default main input of this pipeline is a text file with accession numbers
for the ``reads_download`` component. The ``extra_input`` creates
a new parameter, named ``local`` in this example, that allows us to provide
additional input data to the ``integrity_coverage`` component directly::

    nextflow run pipe.nf --accessions accession_list.txt --local "fastq/*_{1,2}.*"

What will happen in this pipeline, is that the fastq files provided to the
``integrity_coverage`` component will be mixed with the ones provided by the
``reads_download`` component. Therefore, if we provide 10 accessions and 10
fastq samples, we'll end up with 20 samples being processed by the end of the
pieline.

**It is important to note that the extra input parameter expected data
compliant with the input type of the process.** If files other than fastq files
would be provided in the pipeline above, this would result in a pipeline error.

If the ``extra_input`` directive is used on a component that has a different
input type from the first component in the pipeline, it is possible to use
the ``default`` value::

    trimmomatic spades abricate={'extra_input':'default'}

In this case, the input type of the first component if fastq and the input
type of ``abricate`` is fasta. The ``default`` value will make available the
default parameter for fasta raw input, which is ``fasta``::

    nextflow run pipe.nf --fastq "fastq/*_{1,2}.*" --fasta "fasta/*.fasta"

Pipeline file
-------------

Instead of providing the pipeline components via the command line, you can
specify them in a text file::

    # my_pipe.txt
    trimmomatic fastqc spades

And then provide the pipeline file to the ``-t`` parameter::

    assemblerflow build -t my_pipe.txt -o my_pipe.nf

Pipeline files are usually more readable, particularly when they become more
complex. Consider the following example::

    integrity_coverage (
        spades={'memory':'\'50GB\''} |
        skesa={'memory':'\'40GB\'','cpus':'4'} |
        trimmomatic fastqc (
            spades pilon (abricate={'extra_input':'default'} | prokka) |
            skesa pilon (abricate | prokka)
        )
    )

In addition to be more readable, it is also easier to edit, re-use and share.

