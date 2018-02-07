#!/usr/bin/env python3

__version__ = "1.0.0"
__build__ = "22012018"

import os
import shutil
import logging
import argparse
import logging.config

from distutils.dir_util import copy_tree
from os.path import join, dirname

try:
    from generator import HeaderSkeleton as hs
    import generator.Process as pc
except ImportError:
    from assemblerflow.generator import HeaderSkeleton as hs
    import assemblerflow.generator.Process as pc

logger = logging.getLogger("main")


class ProcessError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ChannelError(Exception):
    def __init__(self, p1, p2, t1, t2):
        self.p1 = p1
        self.p2 = p2
        self.t1 = t1
        self.t2 = t2

    def __str__(self):
        return "The output of the '{}' process ({}) cannot link with the " \
               "input of the '{}' process ({}). Please check the order of " \
               "the processes".format(self.p1, self.p2, self.t1, self.t2)


class NextflowGenerator:

    process_map = {
        "integrity_coverage": pc.IntegrityCoverage,
        "seq_typing": pc.SeqTyping,
        "patho_typing": pc.PathoTyping,
        "check_coverage": pc.CheckCoverage,
        "fastqc": pc.FastQC,
        "trimmomatic": pc.Trimmomatic,
        "fastqc_trimmomatic": pc.FastqcTrimmomatic,
        "spades": pc.Spades,
        "process_spades": pc.ProcessSpades,
        "assembly_mapping": pc.AssemblyMapping,
        "pilon": pc.Pilon,
        "mlst": pc.Mlst,
        "abricate": pc.Abricate,
        "prokka": pc.Prokka,
        "chewbbaca": pc.Chewbbaca,
        "status_compiler": pc.StatusCompiler,
        "trace_compiler": pc.TraceCompiler
    }
    """
    dict: Maps the process ids to the corresponding template interface class
    """

    def __init__(self, process_list, nextflow_file):

        # Check if all specified processes are available
        for p in process_list:
            if p not in self.process_map:
                raise ValueError(
                    "The process '{}' is not available".format(p))

        init_process = [pc.Init(template="init")]

        processes = [
            self.process_map[p](template=p) for p in process_list
        ]
        self.processes = init_process + processes
        """
        list: Stores the process interfaces in the specified order
        """

        self.nf_file = nextflow_file
        """
        str: Path to file where the pipeline will be generated
        """

        self.template = ""
        """
        str: String that will harbour the pipeline code
        """

        self.secondary_channels = {}
        """
        dict: Stores secondary channel links
        """

        self.secondary_inputs = {}
        """
        dict: Stores the secondary input channels that may be required by
        some processes. The key is the params variable and the key is the
        channel definition for nextflow::

            {"genomeSize": "IN_genome_size = Channel.value(params.genomeSize)"}

        """

        self.status_channels = []
        """
        list: Stores the status channels from each process
        """

        self._check_pipeline_requirements()

    def _check_pipeline_requirements(self):
        """ Checks for some pipeline requirements before building

        Currently, the only hard requirement is that the pipeline must start
        with the integrity_coverage process, in order to evaluate if the
        input FastQ are corrupt or not.

        Besides this requirements, it checks for the existence the dependencies
        for all processes.
        """

        pipeline_names = [x.template for x in self.processes]

        logger.debug("Checking pipeline requirements for template "
                     "list: {}".format(pipeline_names))

        # Check if the pipeline contains at least one process with raw input
        # type
        raw_processes = [p for p in self.processes if p.input_type == "raw"]
        if not raw_processes:
            raise ProcessError("At least one process with 'raw' input type "
                               "must be specified. Check if the "
                               "pipeline starts with an appropriate starting"
                               " process.")

        logger.debug("Checking for dependencies of templates")

        for p in [i for i in self.processes if i.dependencies]:
            if not set(p.dependencies).issubset(set(pipeline_names)):
                raise ProcessError(
                    "Missing dependencies for process {}: {}".format(
                        p.template, p.dependencies))

    def _build_header(self):
        """Adds the header template to the master template string
        """

        logger.debug("Building header")
        self.template += hs.header

    def _set_channels(self):
        """Sets the main channels for the pipeline

        The setup of the main channels follows four main steps for each
        process specified in the :py:attr:`NextflowGenerator.processes`
        attribute:

            - (If not the first process) Checks if the input of the current
            process is compatible with the output of the previous process.
            - Checks if the current process has starts any secondary channels.
            If so, populate the :py:attr:`NextflowGenerator.secondary_channels`
            with the name of the link start, the process class and a list
            to harbour potential receiving ends.
            - Checks if the current process receives from any secondary
            channels. If a corresponding secondary link has been previously
            set, it will populate the
            :py:attr:`NextflowGenerator.secondary_channels` attribute with
            the receiving channels.
            - Sets the main channels by providing the process ID.

        Notes
        -----
        **On the secondary channel setup**: With this approach, there can only
        be one secondary link start for each type of secondary link. For
        instance, If there are two processes that start a secondary channel
        for the ``SIDE_max_len`` channel, only the last one will be recorded,
        and all receiving processes will get the channel from the latest
        process.
        """

        logger.debug("Setting main channels")
        previous_channel = None
        pidx = 0

        for idx, p in enumerate(self.processes):

            # Make sure that the process id starts at 1
            if not p.ignore_pid:
                pidx += 1
            else:
                logger.debug("[{}] Ignoring process id increment".format(
                    p.template
                ))

            if p.ptype == "terminal":
                # Get last main channel
                channel_str = previous_channel._main_out_str
                p._main_in_str = channel_str
                p.link_end = [{"link": channel_str,
                               "alias": channel_str}]

            logger.debug("[{}] Setting main channels for idx '{}'".format(
                p.template, pidx))
            logger.debug("[{}] Expected input type: {}".format(
                p.template, p.input_type))

            if not previous_channel:
                # Set the first output type
                previous_channel = p
            else:
                logger.debug(
                    "[{}] Previous output type for template: {}".format(
                        p.template, previous_channel.output_type))
                # Check if the connecting processes can be linked by their
                # input/output types
                if p.ignore_type:
                    pass
                elif previous_channel.output_type != p.input_type:
                    raise ChannelError(previous_channel.template,
                                       previous_channel.output_type,
                                       p.template,
                                       p.input_type)
                else:
                    previous_channel = p

            logger.debug("[{}] Checking secondary links".format(p.template))

            # Check if the current process has a start of a secondary
            # side channel
            if p.link_start:
                logger.debug("[{}] Found secondary link start: {}".format(
                    p.template, p.link_start))
                for l in p.link_start:
                    self.secondary_channels[l] = {"p": p, "end": []}

            # check if the current process receives a secondary side channel.
            # If so, add to the links list of that side channel
            if p.link_end:
                logger.debug("[{}] Found secondary link end: {}".format(
                    p.template, p.link_end))
                for l in p.link_end:
                    if l["link"] in self.secondary_channels:
                        self.secondary_channels[l["link"]]["end"].append(
                            "{}_{}".format(l["alias"], pidx))

            if p.status_channels:
                logger.debug("[{}] Added status channel(s): {}".format(
                    p.template, p.status_channels))
                self.status_channels.append(p.status_strs)

            if p.secondary_inputs:
                logger.debug("[{}] Found secondary input channel(s): "
                             "{}".format(p.template, p.secondary_inputs))
                for ch in p.secondary_inputs:
                    if ch["params"] not in self.secondary_inputs:
                        logger.debug("[{}] Added channel: {}".format(
                            p.template, ch["channel"]))
                        self.secondary_inputs[ch["params"]] = ch["channel"]

            logger.debug("[{}] Setting main channels with pid '{}'".format(
                p.template, pidx))

            p.set_channels(**{"pid": pidx})

    def _set_secondary_inputs(self):

        # Get init process
        init_process = self.processes[0]
        init_process.set_secondary_inputs(self.secondary_inputs)

    def _set_secondary_channels(self):
        """Sets the secondary channels for the pipeline

        This will iterate over the
        :py:attr:`NextflowGenerator.secondary_channels` dictionary that is
        populated when executing :py:func:`NextflowGenerator._set_channels`
        method.
        """

        logger.debug("Setting secondary channels: {}".format(
            self.secondary_channels))

        for source, vals in self.secondary_channels.items():

            # Ignore status processes
            if vals["p"].ptype == "status":
                logger.debug("Skipping template {} of type {}".format(
                    vals["p"].template, vals["p"].ptype))
                continue

            # Skip if there are no receiving ends for this secondary channel
            if not vals["end"]:
                logger.debug("[{}] No secondary links to setup".format(
                    vals["p"].template))
                continue

            logger.debug("[{}] Setting secondary links for "
                         "source {}: {}".format(vals["p"].template,
                                                source,
                                                vals["end"]))

            vals["p"].set_secondary_channel(source, vals["end"])

    def _set_status_channels(self):
        """Compiles all status channels for the status compiler process
        """

        # Compile status channels from pipeline process
        status_channels = []
        for p in [p for p in self.processes if p.ptype != "status"]:
            status_channels.extend(p.status_strs)

        logger.debug("Setting status channels: {}".format(status_channels))

        # Check for duplicate channels. Raise exception if found.
        if len(status_channels) != len(set(status_channels)):
            raise ProcessError(
                "Duplicate status channels detected. Please ensure that "
                "the 'status_channels' attributes of each process are "
                "unique. Here are the status channels:\n\n{}".format(
                    ", ".join(status_channels)
                ))

        for p in self.processes:
            if p.ptype == "status":
                p.set_status_channels(status_channels)

    def build(self):
        """Main pipeline builder

        This method is responsible for building the
        :py:attr:`NextflowGenerator.template` attribute that will contain
        the nextflow code of the pipeline.

        First it builds the header, then sets the main channels, the
        secondary channels and finally the status channels. When the pipeline
        is built, is writes the code to a nextflow file.
        """

        # Generate regular nextflow header that sets up the shebang, imports
        # and all possible initial channels
        self._build_header()

        self._set_channels()

        self._set_secondary_inputs()

        self._set_secondary_channels()

        self._set_status_channels()

        for p in self.processes:
            self.template += p.template_str

        with open(self.nf_file, "w") as fh:
            fh.write(self.template)


def get_args():

    parser = argparse.ArgumentParser(
        description="Nextflow pipeline generator")

    parser.add_argument("-t", "--tasks", nargs="+", dest="tasks",
                        help="Space separated tasks of the pipeline")
    parser.add_argument("-o", dest="output_nf",
                        help="Name of the pipeline file")
    parser.add_argument("--include-templates", dest="include_templates",
                        action="store_const", const=True,
                        help="This will copy the necessary templates and lib"
                             " files to the directory where the nextflow"
                             " pipeline will be generated")
    parser.add_argument("--debug", dest="debug", action="store_const",
                        const=True, help="Set log to debug mode")

    args = parser.parse_args()

    return args


def copy_project(path):
    """

    Parameters
    ----------
    path

    Returns
    -------

    """

    # Get nextflow repo directory
    repo_dir = dirname(os.path.abspath(__file__))

    # Get target directory
    target_dir = dirname(path)

    # Copy templates
    copy_tree(join(repo_dir, "templates"), join(target_dir, "templates"))

    # Copy Helper scripts
    copy_tree(join(repo_dir, "lib"), join(target_dir, "lib"))

    # Copy bin scripts
    copy_tree(join(repo_dir, "bin"), join(target_dir, "bin"))

    # Copy default config file
    shutil.copy(join(repo_dir, "nextflow.config"),
                join(target_dir, "nextflow.config"))


def run(args):

    if args.debug:
        logger.setLevel(logging.DEBUG)

        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        # create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # add formatter to ch
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    # nfg = NextflowGenerator(args.tasks, args.output_nf)
    nfg = NextflowGenerator(process_list=args.tasks,
                            nextflow_file=args.output_nf)
    # nfg = NextflowGenerator(pipeline, "/home/diogosilva/teste/teste.nf")

    nfg.build()

    if args.include_templates:
        copy_project(args.output_nf)


def main():

    args = get_args()
    run(args)


if __name__ == '__main__':

    main()
