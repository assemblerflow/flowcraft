Template creation guidelines
============================

Though none of these guidelines are mandatory nor required, their usage is
highly recommended for several reasons:

- Consistency in the outputs of the templates throughout the pipeline,
  particularly the status and report dotfiles (see [dotfiles section]);
- Debugging purposes;
- Versioning;
- Proper documentation of the template scripts.

Preface header
--------------

After the script shebang, a header with a brief description of the purpose and
expected inputs and outputs should be provided. A complete example of such
description can be viewed in :mod:`assemblerflow.templates.integrity_coverage`.

Purpose
^^^^^^^

Purpose section contains a brief description of the script's objective. E.g.::

    Purpose
    -------

    This module is intended parse the results of FastQC for paired end FastQ \
    samples.

Expected input
^^^^^^^^^^^^^^

Expected input section contains a description of the variables that are
provided to the main function of the template script. These variables are
defined in the input channels of the process in which the template is supposed
to be executed. E.g.::

    Expected input
    --------------

    The following variables are expected whether using NextFlow or the
    :py:func:`main` executor.

    - ``mash_output`` : String with the name of the mash screen output file.
        - e.g.: ``'sortedMashScreenResults_SampleA.txt'``

This means that the process that will execute this channel will have the input
defined as::

    input:
    file(mash_output) from <channel>

Generated output
^^^^^^^^^^^^^^^^

Generated output section contains a description of the output files that the
template script is intended to generated. E.g.::

    Generated output
    ----------------

    The generated output are output files that contain an object, usually a string.

    - ``fastqc_health`` : Stores the health check for the current sample. If it
        passes all checks, it contains only the string 'pass'. Otherwise, contains
        the summary categories and their respective results


Versioning and logging
----------------------

Since assemblerflow has a specific ``logger`` and version system, a
requirement should be imported from `templates.utils
<https://github.com/ODiogoSilva/templates/tree/master/utils>`_::

    # the module that imports the logger and the decorator class for versioning
    # of the script itself and other software used in the script
    from utils.assemblerflow_base import get_logger, MainWrapper



Logger
^^^^^^

A `logger` function is also required to add logs to the script.

First, the logger must be called, for example, after the **imports** as follows::

    logger = get_logger(__file__)

Then, it may be used at will, using the default `logging levels
<https://docs.python.org/3.6/library/logging.html#levels>`_ . E.g.::

    logger.error("Module exited unexpectedly with error:\\n{}".format(
                traceback.format_exc()))

MainWrapper decorator
^^^^^^^^^^^^^^^^^^^^^

This class decorator allows the program to fetch information on the script version,
build and template name. For example::

    # This can also be declared after the imports
    __version__ = "1.0.0"
    __build__ = "15012018"
    __template__ = "process_abricate-nf"

The ``MainWrapper`` decorator should be added to the main function of the script.
E.g.::

    @MainWrapper
    def main():
        #some awesome code
        ...

Besides searching for the script's version, build and template name this decorator
will also search for a specific set of functions that start with the
substring ``__set_version``. For example::

    def __set_version_fastqc():

        try:

        cli = ["fastqc", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        version = stdout.strip().split()[1][1:].decode("utf8")

        except Exception as e:
            logger.debug(e)
            version = "undefined"

        # Note that it returns a dictionary that will then be written to the .versions
        # dotfile
        return {
            "program": "FastQC",
            "version": version,
            # some programs may also contain build.
        }


Nextflow `.command.sh`
----------------------

When these templates are used with Nextflow `template <https://www.nextflow.io/docs/latest/process.html#template>`_
a ``.command.sh`` file will be generated, allowing to pass arguments between nextflow
 pipeline and python scripts. In this case, it is recommended that
an **if statement** is included to parse the arguments from nextflow to python template.
For example, imagine we have a path to a file name to pass as argument between
nextflow and the required template::

    # code check for nextflow execution
    if __file__.endswith(".command.sh"):
        FILE_NAME = '$Nextflow_file_name'
        # logger output can also be included here, for example:
        logger.debug("Running {} with parameters:".format(
            os.path.basename(__file__)))
        logger.debug("FILE_NAME: {}".format(FILE_NAME))

Then, we could use this variable as the argument of a function, such as::

    def main(FILE_NAME):
        #some awesome code
        ...


This way, we can use this function with nextflow arguments or without them.

Use numpy docstrings
--------------------

``Assemblerflow`` uses numpy docstrings to document code.
Use
`this link <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_
for an example.