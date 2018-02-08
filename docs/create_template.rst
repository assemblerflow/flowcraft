Template creation guidelines
============================

Preface header
--------------

After the script shebang, a header with a brief description of the purpose and
expected inputs and outputs should be provided:

1. Purpose section contains a brief description of the script's objective. E.g.::

    Purpose
    -------

    This module is intended parse the results of FastQC for paired end FastQ \
    samples.

2. Expected input section contains a description of the variables that are
provided to the main function of the template script. These variables are
defined in the input channels of the process in which the template is supposed
to be executed. E.g.::

    Expected input
    --------------

    The following variables are expected whether using NextFlow or the
    :py:func:`main` executor.

    - ``mash_output`` : String with the name of the mash screen output file.
        - e.g.: ``'sortedMashScreenResults_SampleA.txt'``

3. Generated output section contains a description of the output files that the
template script is intended to generated. E.g.::

    Generated output
    ----------------

    The generated output are output files that contain an object, usually a string.

    - ``fastqc_health`` : Stores the health check for the current sample. If it
        passes all checks, it contains only the string 'pass'. Otherwise, contains
        the summary categories and their respective results


Mandatory requirements
----------------------

Since assemblerflow has a specific `logger`, a set of requirements are required
so that the logger can properly work::

    # standard python packages
    import os
    import sys
    import traceback

    # try/except used to search for the utils.assemblerflow in path if used
    # within a docker image
    try:
        sys.path.append(os.environ["ASSEMBLERFLOW_UTILS"])
    except KeyError:
        pass
    # import the logger it self
    from utils.assemblerflow_base import get_logger, _log_error

Then, the logger must be called as::

    logger = get_logger(__file__)

And finally you may use the logger as you want, using the default `logging levels
<https://docs.python.org/3.6/library/logging.html#levels>`_ . E.g.::

    logger.error("Module exited unexpectedly with error:\\n{}".format(
                traceback.format_exc()))


Checks for versions (build_versions)
------------------------------------

A `build_versions` function which has the versions of the script and programs
used by the template script and that can be logged using the `logger` generated
in `Mandatory requirements`_. E.g.::

    def build_versions():
        logger.debug("Checking module versions")

        ver = [{
            "program": __template__,
            "version": __version__,
            "build": __build__
        }]
        logger.debug("Versions list set to: {}".format(ver))

        with open(".versions", "w") as fh:
            fh.write(json.dumps(ver, separators=(",", ":")))

Other programs versions can also be added to the `ver` variable, by adding a
function that obtains this information from shell using `subprocess`. E.g.::

    def get_abricate_version():

        try:

            # Get abricate version
            cli = ["abricate", "--version"]
            p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
            stdout, _ = p.communicate()

            version = stdout.strip().split()[-1].decode("utf8")

        except Exception as e:
            logger.debug(e)
            version = "undefined"

        try:

            # Get abricate database versions
            cli = ["abricate", "--list"]
            p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
            dbout, _ = p.communicate()

            databases = [[u.decode("utf8") for u in i.strip().split()]
                         for i in dbout.splitlines()][1:]

        except Exception as e:
            logger.debug(e)
            databases = "undefined"

        return {
            "program": "abricate",
            "version": version,
            "databases": databases
        }

Try/except block in main execution
----------------------------------

A try except block in main execution is required so that an error can be raised
if something goes really wrong with the template script execution. E.g.::

    try:
        build_versions()
        main(FASTQ_ID, RESULT_P1, RESULT_P2, OPTS)
    except:
        logger.error("Module exited unexpectedly with error:\\n{}".format(
            traceback.format_exc()))
        _log_error()

Dotfiles
--------

.status
^^^^^^^

.warning
^^^^^^^^

.fail
^^^^^

.report.json
^^^^^^^^^^^^

.versions
^^^^^^^^^
