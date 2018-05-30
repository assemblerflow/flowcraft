"""

"""

import os
import sys
import json
import logging
import traceback

from time import gmtime, strftime


def get_logger(filepath, level=logging.DEBUG):
    # create logger
    logger = logging.getLogger(os.path.basename(filepath))
    logger.setLevel(level)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(level)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)

    return logger


def log_error():
    """Nextflow specific function that logs an error upon unexpected failing
    """

    with open(".status", "w") as status_fh:
        status_fh.write("error")


class MainWrapper:

    def __init__(self, f):

        self.f = f
        self.context = self.f.__globals__
        self.logger = self.context.get("logger", None)

    def __call__(self, *args, **kwargs):

        self.logger.debug("Starting template at {}".format(
            strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        self.logger.debug("Working directory: {}".format(os.getcwd()))

        try:
            self.build_versions()
            self.f(*args, **kwargs)
        except SystemExit as e:
            sys.exit(e)
        except:
            if self.logger:
                self.logger.error("Module exited unexpectedly with error:"
                                  "\\n{}".format(traceback.format_exc()))
            log_error()

        self.logger.debug("Finished template at {}".format(
            strftime("%Y-%m-%d %H:%M:%S", gmtime())))

    def build_versions(self):
        """Writes versions JSON for a template file

        This method creates the JSON file ``.versions`` based on the metadata
        and specific functions that are present in a given template script.

        It starts by fetching the template metadata, which can be specified
        via the ``__version__``, ``__template__`` and ``__build__``
        attributes. If all of these attributes exist, it starts to populate
        a JSON/dict array (Note that the absence of any one of them will
        prevent the version from being written).

        Then, it will search the
        template scope for functions that start with the substring
        ``__set_version`` (For example ``def __set_version_fastqc()`).
        These functions should gather the version of
        an arbitrary program and return a JSON/dict object with the following
        information::

            {
                "program": <program_name>,
                "version": <version>
                "build": <build>
            }

        This JSON/dict object is then written in the ``.versions`` file.
        """

        version_storage = []

        template_version = self.context.get("__version__", None)
        template_program = self.context.get("__template__", None)
        template_build = self.context.get("__build__", None)

        if template_version and template_program and template_build:
            if self.logger:
                self.logger.debug("Adding template version: {}; {}; "
                                  "{}".format(template_program,
                                              template_version,
                                              template_build))
            version_storage.append({
                "program": template_program,
                "version": template_version,
                "build": template_build
            })

        for var, obj in self.context.items():
            if var.startswith("__get_version"):
                ver = obj()
                version_storage.append(ver)
                if self.logger:
                    self.logger.debug("Found additional software version"
                                      "{}".format(ver))

        with open(".versions", "w") as fh:
            fh.write(json.dumps(version_storage, separators=(",", ":")))

