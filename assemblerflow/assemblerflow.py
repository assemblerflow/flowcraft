#!/usr/bin/env python3

__version__ = "1.0.0"
__build__ = "22012018"

import os
import sys
import shutil
import logging
import argparse
import logging.config

from distutils.dir_util import copy_tree
from collections import defaultdict
from os.path import join, dirname

try:
    from generator import HeaderSkeleton as hs
    from generator.pipeline_parser import parse_pipeline, SanityError
    from generator.process_details import proc_collector, colored_print
    import generator.Process as pc
except ImportError:
    from assemblerflow.generator import HeaderSkeleton as hs
    from assemblerflow.generator.pipeline_parser import parse_pipeline, \
        SanityError
    from assemblerflow.generator.process_details import proc_collector, \
        colored_print
    import assemblerflow.generator.Process as pc

logger = logging.getLogger("main")

process_map = {
        "integrity_coverage": pc.IntegrityCoverage,
        "seq_typing": pc.SeqTyping,
        "patho_typing": pc.PathoTyping,
        "check_coverage": pc.CheckCoverage,
        "fastqc": pc.FastQC,
        "trimmomatic": pc.Trimmomatic,
        "fastqc_trimmomatic": pc.FastqcTrimmomatic,
        "skesa": pc.Skesa,
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

    def __init__(self, process_list, nextflow_file):

        # Check if all specified processes are available
        for p in process_list:
            pname = p["input"]["process"]
            if pname not in process_map and pname != "__init__":
                logger.error("The process '{}' is not available".format(pname))
                sys.exit(1)

        self.processes = []

        self.processes = [pc.Init(template="init")]
        """
        list: Stores the process interfaces in the specified order
        """

        self._fork_tree = defaultdict(list)

        self._build_connections(process_list)

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

        self.main_raw_inputs = {}
        """
        list: Stores the main raw inputs from the user parameters into the
        first process(es). 
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

        # self._check_pipeline_requirements()

    def _build_connections(self, process_list):
        """

        Returns
        -------

        """

        logger.debug("=============================")
        logger.debug("Building pipeline connections")
        logger.debug("=============================")

        for p, con in enumerate(process_list):

            logger.debug("Processing connection '{}': {}".format(p, con))

            # Get lanes
            in_lane = con["input"]["lane"]
            out_lane = con["output"]["lane"]
            logger.debug("[{}] Input lane: {}".format(p, in_lane))
            logger.debug("[{}] Output lane: {}".format(p, out_lane))

            p_in_name = con["input"]["process"]
            logger.debug("[{}] Input channel: {}".format(p, p_in_name))
            p_out_name = con["output"]["process"]
            logger.debug("[{}] Output channel: {}".format(p, p_out_name))

            # Instance output process
            out_process = process_map[p_out_name](template=p_out_name)
            input_suf = "{}_{}".format(in_lane, p)
            output_suf = "{}_{}".format(out_lane, p)
            logger.debug("[{}] Setting main channels with input suffix '{}'"
                         " and output suffix '{}'".format(p, input_suf,
                                                          output_suf))
            out_process.set_main_channel_names(input_suf, output_suf, out_lane)

            # Instance input process, if it exists. In case of init, the
            # output process forks from the raw input user data
            if p_in_name != "__init__":
                in_process = process_map[p_in_name](template=p_in_name)
                # Test if two processes can be connected by input/output types
                logger.debug("[{}] Testing connection between input and "
                             "output processes".format(p))
                self._test_connection(in_process, out_process)
                out_process.parent_lane = in_lane
            else:
                out_process.parent_lane = None

            # If the current connection is a fork, add it to the fork tree
            if in_lane != out_lane:
                logger.debug("[{}] Connection is a fork. Adding lanes to "
                             "fork list".format(p))
                self._fork_tree[in_lane].append(out_lane)
                # Update main output fork of parent process
                try:
                    parent_fork = [x for x in self.processes
                                   if x.lane == in_lane and
                                   x.template == p_in_name][0]
                    logger.debug("[{}] Updating main forks of parent fork "
                                 "'{}' with '{}'".format(
                                    p, parent_fork, out_process.input_channel))
                    parent_fork.update_main_forks(out_process.input_channel)
                except IndexError:
                    pass
            else:
                parent_process = self.processes[-1]
                if parent_process.output_channel:
                    logger.debug(
                        "[{}] Updating input channel of output process"
                        " with '{}'".format(
                            p, parent_process.output_channel))
                    out_process.input_channel = parent_process.output_channel

            self.processes.append(out_process)

    @staticmethod
    def _test_connection(parent_process, child_process):
        """Tests if two processes can be connected by input/output type

        Parameters
        ----------
        parent_process : assemblerflow.Process.Process
            Process that will be sending output.
        child_process : assemblerflow.Process.Process
            Process that will receive output.

        """

        # If any of the processes has an ignore type attribute set to True,
        # don't perform the check
        if parent_process.ignore_type or child_process.ignore_type:
            return

        if parent_process.output_type != child_process.input_type:
            logger.error(
                "The output of the '{}' process ({}) cannot link with the "
                "input of the '{}' process ({}). Please check the order of "
                "the processes".format(parent_process.template,
                                       parent_process.output_type,
                                       child_process.template,
                                       child_process.input_type))
            sys.exit(1)

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

    def _update_raw_input(self, p):
        """Given a process, this method updates the
        :attr:`~Process.main_raw_inputs` attribute with the corresponding
        raw input channel of that process

        Parameters
        ----------
        p : assemblerflow.Process.Process
        """

        logger.debug("[{}] Setting raw input channel".format(p.template))
        raw_in = p.get_user_channel()

        if p.input_type in self.main_raw_inputs:
            self.main_raw_inputs[p.input_type]["raw_forks"].append(
                raw_in["input_channel"])
        else:
            self.main_raw_inputs[p.input_type] = {
                "channel": raw_in["channel"],
                "channel_str": raw_in["channel_str"],
                "raw_forks": [raw_in["input_channel"]]
            }

    def _update_secondary_inputs(self, p):
        """Given a process, this method updates the
        :attr:`~Process.secondary_inputs` attribute with the corresponding
        secondary inputs of that channel.

        Parameters
        ----------
        p : assemblerflow.Process.Process
        """

        logger.debug("[{}] Checking secondary links".format(p.template))
        if p.secondary_inputs:
            logger.debug("[{}] Found secondary input channel(s): "
                         "{}".format(p.template, p.secondary_inputs))
            for ch in p.secondary_inputs:
                if ch["params"] not in self.secondary_inputs:
                    logger.debug("[{}] Added channel: {}".format(
                        p.template, ch["channel"]))
                    self.secondary_inputs[ch["params"]] = ch["channel"]

    def _get_fork_tree(self, p):
        """

        Parameters
        ----------
        p

        Returns
        -------
        """

        lane = p.lane
        parent_lanes = [lane]

        while True:
            original_lane = lane
            for fork_in, fork_out in self._fork_tree.items():
                if lane in fork_out:
                    lane = fork_in
                    parent_lanes.append(fork_in)
            if lane == original_lane:
                break

        return parent_lanes

    def _update_secondary_channels(self, p):
        """

        Parameters
        ----------
        p : assemblerflow.Process.Process

        """

        # Check if the current process has a start of a secondary
        # side channel
        if p.link_start:
            logger.debug("[{}] Found secondary link start: {}".format(
                p.template, p.link_start))
            for l in p.link_start:
                self.secondary_channels[l] = {p.lane: {"p": p, "end": []}}

        # check if the current process receives a secondary side channel.
        # If so, add to the links list of that side channel
        if p.link_end:
            logger.debug("[{}] Found secondary link end: {}".format(
                p.template, p.link_end))
            for l in p.link_end:

                if l["link"] not in self.secondary_channels:
                    continue

                parent_forks = self._get_fork_tree(p)
                for lane in parent_forks:
                    if lane in self.secondary_channels[l["link"]]:
                        self.secondary_channels[
                            l["link"]][lane]["end"].append("{}".format(
                                "{}_{}".format(l["alias"], p.pid)))

        logger.debug("[{}] Secondary links updated: {}".format(
            p.template, self.secondary_channels))

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

        for i, p in enumerate(self.processes):

            # Set main channels for the process
            logger.debug("[{}] Setting main channels with pid: {}".format(
                p.template, i))
            p.set_channels(pid=i)

            # If there is no parent lane, set the raw input channel from user
            if not p.parent_lane and p.input_type:
                self._update_raw_input(p)

            self._update_secondary_inputs(p)

            self._update_secondary_channels(p)

        # for idx, p in enumerate(self.processes):
        #
        #     # Make sure that the process id starts at 1
        #     if not p.ignore_pid:
        #         pidx += 1
        #     else:
        #         logger.debug("[{}] Ignoring process id increment".format(
        #             p.template
        #         ))
        #
        #     if p.ptype == "terminal":
        #         # Get last main channel
        #         channel_str = previous_channel._main_out_str
        #         p._main_in_str = channel_str
        #         p.link_end = [{"link": channel_str,
        #                        "alias": channel_str}]
        #
        #     logger.debug("[{}] Setting main channels for idx '{}'".format(
        #         p.template, pidx))
        #     logger.debug("[{}] Expected input type: {}".format(
        #         p.template, p.input_type))
        #
        #     if not previous_channel:
        #         # Set the first output type
        #         previous_channel = p
        #     else:
        #         logger.debug(
        #             "[{}] Previous output type for template: {}".format(
        #                 p.template, previous_channel.output_type))
        #         # Check if the connecting processes can be linked by their
        #         # input/output types
        #         if p.ignore_type:
        #             pass
        #         elif previous_channel.output_type != p.input_type:
        #             raise ChannelError(previous_channel.template,
        #                                previous_channel.output_type,
        #                                p.template,
        #                                p.input_type)
        #         else:
        #             previous_channel = p
        #
        #     logger.debug("[{}] Checking secondary links".format(p.template))
        #
        #     # Check if the current process has a start of a secondary
        #     # side channel
        #     if p.link_start:
        #         logger.debug("[{}] Found secondary link start: {}".format(
        #             p.template, p.link_start))
        #         for l in p.link_start:
        #             self.secondary_channels[l] = {"p": p, "end": []}
        #
        #     # check if the current process receives a secondary side channel.
        #     # If so, add to the links list of that side channel
        #     if p.link_end:
        #         logger.debug("[{}] Found secondary link end: {}".format(
        #             p.template, p.link_end))
        #         for l in p.link_end:
        #             if l["link"] in self.secondary_channels:
        #                 self.secondary_channels[l["link"]]["end"].append(
        #                     "{}_{}".format(l["alias"], pidx))
        #
        #     if p.status_channels:
        #         logger.debug("[{}] Added status channel(s): {}".format(
        #             p.template, p.status_channels))
        #         self.status_channels.append(p.status_strs)
        #
        #     if p.secondary_inputs:
        #         logger.debug("[{}] Found secondary input channel(s): "
        #                      "{}".format(p.template, p.secondary_inputs))
        #         for ch in p.secondary_inputs:
        #             if ch["params"] not in self.secondary_inputs:
        #                 logger.debug("[{}] Added channel: {}".format(
        #                     p.template, ch["channel"]))
        #                 self.secondary_inputs[ch["params"]] = ch["channel"]
        #
        #     logger.debug("[{}] Setting main channels with pid '{}'".format(
        #         p.template, pidx))
        #
        #     p.set_channels(**{"pid": pidx})

    def _set_secondary_inputs(self):

        # Get init process
        init_process = self.processes[0]
        init_process.set_raw_inputs(self.main_raw_inputs)
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

        for source, lanes in self.secondary_channels.items():

            for lane, vals in lanes.items():

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

        # self._set_status_channels()

        for p in self.processes:
            self.template += p.template_str

        with open(self.nf_file, "w") as fh:
            fh.write(self.template)


def get_args():

    parser = argparse.ArgumentParser(
        description="Nextflow pipeline generator")

    group_lists = parser.add_mutually_exclusive_group()

    parser.add_argument("-t", "--tasks", type=str, dest="tasks",
                        help="Space separated tasks of the pipeline")
    parser.add_argument("-o", dest="output_nf",
                        help="Name of the pipeline file")
    parser.add_argument("--include-templates", dest="include_templates",
                        action="store_const", const=True,
                        help="This will copy the necessary templates and lib"
                             " files to the directory where the nextflow"
                             " pipeline will be generated")
    parser.add_argument("-c", "--check-pipeline", dest="check_only",
                        action="store_const", const=True,
                        help="Check only the validity of the pipeline"
                             "string and exit.")
    group_lists.add_argument("-L", "--detailed-list", action="store_const",
                        dest="detailed_list", const=True,
                        help="Print a detailed description for all the "
                             "currently available processes")
    group_lists.add_argument("-l", "--short-list", action="store_const",
                        dest="short_list", const=True,
                        help="Print a short list of the currently available "
                             "processes")
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

        # create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    else:
        logger.setLevel(logging.INFO)

        # create special formatter for info logs
        formatter = logging.Formatter('%(message)s')

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # add formatter to ch
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    welcome = [
        "========= A S S E M B L E R F L O W =========\n",
        "version: {}".format(__version__),
        "build: {}\n".format(__build__),
        "============================================="
    ]

    logger.info(colored_print("1;32m", "\n".join(welcome)))

    # prints a detailed list of the process class arguments
    if args.detailed_list:
        # list of attributes to be passed to proc_collector
        arguments_list = [
            "input_type",
            "output_type",
            "description",
            "dependencies",
            "conflicts"
        ]

        proc_collector(process_map, arguments_list)
        sys.exit(0)

    # prints a short list with each process and the corresponding description
    if args.short_list:
        arguments_list = [
            "description"
        ]
        proc_collector(process_map, arguments_list)
        sys.exit(0)

    try:
        logger.info(colored_print(
            "1;38m", "Checking pipeline for errors..."
        ))
        pipeline_list = parse_pipeline(args.tasks)
    except SanityError as e:
        logger.error(colored_print("1;31m", e.value))
        sys.exit(1)
    logger.debug("Pipeline successfully parsed: {}".format(pipeline_list))

    # Exit if only the pipeline parser needs to be checked
    if args.check_only:
        sys.exit()

    nfg = NextflowGenerator(process_list=pipeline_list,
                            nextflow_file=args.output_nf)

    logger.info(colored_print(
        "1;38m", "\nBuilding your awesome pipeline..."
    ))

    # building the actual pipeline nf file
    nfg.build()

    # copy template to cwd, to allow for immediate execution
    if args.include_templates:
        copy_project(args.output_nf)


def main():

    args = get_args()
    run(args)


if __name__ == '__main__':

    main()
