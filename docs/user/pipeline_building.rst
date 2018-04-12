Pipeline building
=================

Pipeline forks
--------------

The output of any component in an assemblerflow pipeline can be forked into
two or more components, using the following fork syntax::

    trimmomatic fastqc ( spades | skesa )

In this example, the output of ``fastqc`` will be fork into two new *lanes*,
which will proceed independently from each other. In this syntax, a fork is
triggered by the ``(`` symbol (and the corresponding closing ``)``) and each
lane will be separated by a ``|`` symbol. There is no limitation to the number
of forks or lanes that a pipeline has. For instance, we could add more
components after the ``skesa`` module, including another fork::

    trimmomatic fastqc ( spades | skesa pilon (abricate | prokka | chewbbaca) )

In this example, data will be forked after ``fastqc`` into two new lanes,
processed by ``spades`` and ``skesa``. In the skesa lane, data will continue
to flow into the ``pilon`` component and its output will fork into three new
lanes.

.. warning::
    Pay special attention to the syntax of the pipeline string when using
    forks. However, when unable to parse it, assemblerflow will do its best
    to inform you where the parsing error occurred.

Directives
----------

Several directives with information on cpu usage, RAM, version, etc. can be
specified to each individual component when building the pipeline using the
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
                            spades={'memory':'\"10GB\"'}" -o my_pipeline.nf

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
    the ``"`` as in this example: ``spades={'memory':'\"10GB\"'}"``. The
    reason why these directives are not automatically converted is to allow
    the specification of dynamic computing resources, such as
    ``spades={'memory':'{10.Gb*task.attempt}'}"``

Extra inputs
------------

Pipeline file
-------------

Instead of providing the pipeline components via the command line, you can
specify them in a text file::

    # my_pipe.txt
    trimmomatic fastqc spades

And then provide the pipeline file to the ``-t`` parameter::

    assemblerflow build -t my_pipe.txt -o my_pipe.nf

This is equivalent to the first build example. Pipeline files are usually more
readable

Pipeline files
==============

