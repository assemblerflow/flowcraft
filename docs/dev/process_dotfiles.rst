.. _dotfiles:

Dotfiles
========

Several dotfiles (files prefixed by a single ``.``, as in ``.status``) are
created at the beginning of every nextflow process that has the following
placeholder (see :ref:`create-process`): ::

    process myProcess {
        {% include "post.txt" ignore missing %}
        (...)
    }

The actual script that creates the dotfiles is found in
``flowcraft/bin``, is called ``set_dotfiles.sh`` and executes the
following command::

    touch .status .warning .fail .report.json .versions

Status
------

The ``.status`` file simply stores a string with the run status of the process.
The supported status are:

- ``pass``: The process finished successfully
- ``fail``: The process ran without unexpected issues but failed due to some
  quality control check
- ``error``: The process exited with an unexpected error.

Warning
-------

The ``.warning`` file stores any warnings that may occur during the execution
of the process. There is no particular format for the warning messages other
than that each individual warning should be in a separate line.

Fail
----

The ``.fail`` file stores any fail messages that may occur during the
execution of the process. When this occurs, the ``.status`` channel must have
the ``fail`` string as well. As in the warning dotfile, there is no
particular format for the fail message.

.. _report-json:

Report JSON
-----------

.. important::
    The general specification of the report JSON changed in version 1.2.2.
    See the `issue tracker <https://github.com/assemblerflow/flowcraft/issues/96>`_
    for details.

The ``.report.json`` file stores any information from a given process that is
deemed worthy of being reported and displayed at the end of the pipeline.
Any information can be stored in this file, as long as it is in JSON format,
but there are a couple of recommendations that are necessary to follow
for them to be processed by a reporting web app (Currently hosted at
`flowcraft-webapp <https://github.com/assemblerflow/flowcraft-webapp>`_). However, if
data processing will be performed with custom scripts, feel free to specify
your own format.

Information for tables
^^^^^^^^^^^^^^^^^^^^^^

Information meant to be displayed in tables should be in the following
format::

    json_dic = {
        "tableRow": [{
            "sample": "A",
            "data": [{
                "header": "Raw BP",
                "value": 123,
                "table": "qc"
            }, {
                "header": "Coverage",
                "value": 32,
                "table": "qc"
            }]
        }, {
            "sample": "B",
            "data": [{
                "header": "Coverage",
                "value": 35,
                "table": "qc"
            }]
        }]
    }

This provides table information for multiple samples in the same process. In
this case, data for two samples is provided. For each sample, values for
one or more headers can be provided. For instance, this report provides
information about the **Raw BP** and **Coverage** for sample **A** and this
information should go to the **qc** table. If any other information is relevant
to build the table, feel free to add more elements to the JSON.

Information for plots
^^^^^^^^^^^^^^^^^^^^^

Information meant to be displayed in plots should be in the following format::

    json_dic = {
        "plotData": [{
            "sample": "strainA",
            "data": {
                "sparkline": 23123,
                "otherplot": [1,2,3]
             }
        }],
    }

As in the table JSON, *plotData* should be an array with an entry for each
sample. The data for each sample should be another JSON where the keys are
the *plot signatures*, so that we know to which plot the data belongs. The
corresponding values are whatever data object you need.

Other information
^^^^^^^^^^^^^^^^^

Other than tables and plots, which have a somewhat predefined format, there
is not particular format for other information. They will simply store the
data of interest to report and it will be the job of a downstream report app
to process that data into an actual visual report.

.. _versions:

Versions
--------

The ``.version`` dotfile should contain a list of JSON objects with the
version information of the programs used in any given process. There are
only two required key:value pairs:

- ``program``: String with the name of the software/script/template
- ``version``: String with the version of said software.

As an example::

    version = {
        "program": "abricate"
        "version": "0.3.7"
    }

Key:value pairs with other metadata can be included at will for downstream
processing.