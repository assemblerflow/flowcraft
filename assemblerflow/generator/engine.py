import sys
import json
import jinja2
import logging

from collections import defaultdict
from os.path import dirname, join, abspath, split


logger = logging.getLogger("main.{}".format(__name__))

try:
    import generator.process as pc
    import generator.error_handling as eh
    from generator import header_skeleton as hs
    from generator import footer_skeleton as fs
    from generator.process_details import colored_print
except ImportError as e:
    import assemblerflow.generator.process as pc
    import assemblerflow.generator.error_handling as eh
    from assemblerflow.generator import header_skeleton as hs
    from assemblerflow.generator import footer_skeleton as fs
    from assemblerflow.generator.process_details import colored_print


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
        # "status_compiler": pc.StatusCompiler,
        # "trace_compiler": pc.TraceCompiler
}
"""
dict: Maps the process ids to the corresponding template interface class wit
the format::
    
    {
        "<template_string>": pc.TemplateClass
    }
"""


class NextflowGenerator:

    def __init__(self, process_connections, nextflow_file):

        self.processes = []

        # Create the processes attribute with the first special init process.
        # This process will handle the forks of the raw input channels and
        # secondary inputs
        self.processes = [pc.Init(template="init")]
        """
        list: Stores the process interfaces in the specified order
        """

        self._fork_tree = defaultdict(list)
        """
        dict: A dictionary with the fork tree of the pipeline, which consists
        on the the paths of each lane. For instance, a single fork with two
        sinks is represented as: {1: [2,3]}. Subsequent forks are then added
        sequentially: {1:[2,3], 2:[3,4,5]}. This allows the path upstream
        of a process in a given lane to be traversed until the start of the
        pipeline. 
        """

        # Builds the connections in the processes, which parses the
        # process_connections dictionary into the self.processes attribute
        # list.
        self._build_connections(process_connections)

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

        self.skip_class = [pc.Status]
        """
        list: Stores the Process classes that should be skipped when iterating
        over the :attr:`~NextflowGenerator.processes` list.
        """

        self.resources = ""
        """
        str: Stores the resource directives string for each nextflow process. 
        """

        self.containers = ""
        """
        str: Stores the container directives string for each nextflow process.
        """

    @staticmethod
    def _parse_process_name(name_str):
        """Parses the process string and returns the process name and its
        directives

        Process strings my contain directive information with the following
        syntax::

            proc_name={'directive':'val'}

        This method parses this string and returns the process name as a
        string and the directives information as a dictionary.

        Parameters
        ----------
        name_str : str
            Raw string with process name and, potentially, directive
            information

        Returns
        -------
        str
            Process name
        dict or None
            Process directives
        """

        directives = None

        fields = name_str.split("=")
        process_name = fields[0]

        if len(fields) == 2:
            _directives = fields[1].replace("'", '"')
            try:
                directives = json.loads(_directives)
            except json.decoder.JSONDecodeError:
                raise eh.ProcessError(
                    "Could not parse directives for process '{}'. The raw"
                    " string is: {}\n"
                    "Possible causes include:\n"
                    "\t1. Spaces inside directives\n"
                    "\t2. Missing '=' symbol before directives\n"
                    "\t3. Missing quotes (' or \") around directives\n"
                    "A valid example: process_name={{'cpus':'2'}}".format(
                        process_name, name_str))

        return process_name, directives

    def _build_connections(self, process_list):
        """Parses the process connections dictionaries into a process list

        This method is called upon instantiation of the NextflowGenerator
        class. Essentially, it sets the main input/output channel names of the
        processes so that they can be linked correctly.

        If a connection between two consecutive process is not possible due
        to a mismatch in the input/output types, it exits with an error.

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

            # Get process names
            try:
                _p_in_name = con["input"]["process"]
                p_in_name, _ = self._parse_process_name(_p_in_name)
                logger.debug("[{}] Input channel: {}".format(p, p_in_name))
                _p_out_name = con["output"]["process"]
                p_out_name, out_directives = self._parse_process_name(
                    _p_out_name)
                logger.debug("[{}] Output channel: {}".format(p, p_out_name))
            # Exception is triggered when the process name/directives cannot
            # be processes.
            except eh.ProcessError as ex:
                logger.error(colored_print(ex.value, "red_bold"))
                sys.exit(1)

            # Instance output process
            if p_out_name not in process_map:
                logger.error(colored_print(
                    "\nThe process '{}' is not available".format(p_out_name),
                    "red_bold"))
                sys.exit(1)

            out_process = process_map[p_out_name](template=p_out_name)

            # Update directives, if provided
            if out_directives:
                out_process.update_directives(out_directives)

            # Set suffix strings for main input/output channels
            input_suf = "{}_{}".format(in_lane, p)
            output_suf = "{}_{}".format(out_lane, p)
            logger.debug("[{}] Setting main channels with input suffix '{}'"
                         " and output suffix '{}'".format(
                            p, input_suf, output_suf))
            out_process.set_main_channel_names(input_suf, output_suf, out_lane)

            # Instance input process, if it exists. In case of init, the
            # output process forks from the raw input user data
            if p_in_name != "__init__":
                # Create instance of input process
                in_process = process_map[p_in_name](template=p_in_name)
                # Test if two processes can be connected by input/output types
                logger.debug("[{}] Testing connection between input and "
                             "output processes".format(p))
                self._test_connection(in_process, out_process)
                out_process.parent_lane = in_lane
            else:
                # When the input process is __init__, set the parent_lane
                # to None. This will tell the engine that the main input
                # channel the output process of this connection will received
                # from the raw user input.
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

                # Get parent process, naive version
                parent_process = self.processes[-1]

                # Check if the last process' lane matches the lane of the
                # current output process. If not, get the last process
                # in the same lane
                if parent_process.lane and parent_process.lane != out_lane:
                    parent_process = [x for x in self.processes[::-1]
                                      if x.lane == out_lane][0]

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

    def _build_header(self):
        """Adds the header template to the master template string
        """

        logger.debug("===============")
        logger.debug("Building header")
        logger.debug("===============")
        self.template += hs.header

    def _build_footer(self):
        """Adds the footer template to the master template string"""

        logger.debug("===============")
        logger.debug("Building header")
        logger.debug("===============")
        self.template += fs.footer

    def _update_raw_input(self, p, sink_channel=None, input_type=None):
        """Given a process, this method updates the
        :attr:`~Process.main_raw_inputs` attribute with the corresponding
        raw input channel of that process. The input channel and input type
        can be overridden if the `input_channel` and `input_type` arguments
        are provided.

        Parameters
        ----------
        p : assemblerflow.Process.Process
            Process instance whose raw input will be modified
        sink_channel: str
            Sets the channel where the raw input will fork into. It overrides
            the process's `input_channel` attribute.
        input_type: str
            Sets the type of the raw input. It overrides the process's
            `input_type` attribute.
        """

        process_input = input_type if input_type else p.input_type
        process_channel = sink_channel if sink_channel else p.input_channel

        logger.debug("[{}] Setting raw input channel "
                     "with input type '{}'".format(p.template, process_input))
        # Get the dictionary with the raw forking information for the
        # provided input
        raw_in = p.get_user_channel(process_channel, process_input)
        logger.debug("[{}] Fetched process raw user: {}".format(p.template,
                                                                raw_in))

        if process_input in self.main_raw_inputs:
            self.main_raw_inputs[process_input]["raw_forks"].append(
                raw_in["input_channel"])
        else:
            self.main_raw_inputs[process_input] = {
                "channel": raw_in["channel"],
                "channel_str": raw_in["channel_str"],
                "raw_forks": [raw_in["input_channel"]]
            }
        logger.debug("[{}] Updated main raw inputs: {}".format(
            p.template, self.main_raw_inputs))

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

    def _set_implicit_link(self, p, link):
        """

        Parameters
        ----------
        p
        link

        Returns
        -------

        """

        output_type = link["link"].lstrip("_")
        parent_forks = self._get_fork_tree(p)
        fork_sink = "{}_{}".format(link["alias"], p.pid)

        for proc in self.processes[::-1]:
            if proc.lane not in parent_forks:
                continue
            if proc.output_type == output_type:
                proc.update_main_forks(fork_sink)
                logger.debug("[{}] Found special implicit link '{}' with "
                             "output type '{}'. Linked '{}' with process "
                             "{}".format(
                                     p.template, link["link"], output_type,
                                     link["alias"], proc))
                return

        self._update_raw_input(p, fork_sink, output_type)

    def _update_secondary_channels(self, p):
        """Given a process, this method updates the
        :attr:`~Process.secondary_channels` attribute with the corresponding
        secondary inputs of that channel.

        The rationale of the secondary channels is the following:

            - Start storing any secondary emitting channels, by checking the
              `link_start` list attribute of each process. If there are
              channel names in the link start, it adds to the secondary
              channels dictionary.
            - Check for secondary receiving channels, by checking the
              `link_end` list attribute. If the link name starts with a
              `__` signature, it will created an implicit link with the last
              process with an output type after the signature. Otherwise,
              it will check is a corresponding link start already exists in
              the at least one process upstream of the pipeline and if so,
              it will update the ``secondary_channels`` attribute with the
              new link.

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
                # If there are multiple link starts in the same lane, the
                # last one is the only one saved.
                if l in self.secondary_channels:
                    self.secondary_channels[l][p.lane] = {"p": p, "end": []}
                else:
                    self.secondary_channels[l] = {p.lane: {"p": p, "end": []}}

        # check if the current process receives a secondary side channel.
        # If so, add to the links list of that side channel
        if p.link_end:
            logger.debug("[{}] Found secondary link end: {}".format(
                p.template, p.link_end))
            for l in p.link_end:

                # Get list of lanes from the parent forks.
                parent_forks = self._get_fork_tree(p)

                # Parse special case where the secondary channel links with
                # the main output of the specified type
                if l["link"].startswith("__"):
                    self._set_implicit_link(p, l)
                    continue

                # Skip if there is no match for the current link in the
                # secondary channels
                if l["link"] not in self.secondary_channels:
                    continue

                for lane in parent_forks:
                    if lane in self.secondary_channels[l["link"]]:
                        self.secondary_channels[
                            l["link"]][lane]["end"].append("{}".format(
                                "{}_{}".format(l["alias"], p.pid)))

        logger.debug("[{}] Secondary links updated: {}".format(
            p.template, self.secondary_channels))

    def _set_channels(self):
        """Sets the main channels for the pipeline

        This method will parse de the :attr:`~Process.processes` attribute
        and perform the following tasks for each process:

            - Sets the input/output channels and main input forks and adds
              them to the process's
              :attr:`assemblerflow.process.Process._context`
              attribute (See
              :func:`~NextflowGenerator.set_channels`).
            - Automatically updates the main input channel of the first
              process of each lane so that they fork from the user provide
              parameters (See
              :func:`~NextflowGenerator._update_raw_input`).
            - Check for the presence of secondary inputs and adds them to the
              :attr:`~NextflowGenerator.secondary_inputs` attribute.
            - Check for the presence of secondary channels and adds them to the
              :attr:`~NextflowGenerator.secondary_channels` attribute.

        Notes
        -----
        **On the secondary channel setup**: With this approach, there can only
        be one secondary link start for each type of secondary link. For
        instance, If there are two processes that start a secondary channel
        for the ``SIDE_max_len`` channel, only the last one will be recorded,
        and all receiving processes will get the channel from the latest
        process. Secondary channels can only link if the source process if
        downstream of the sink process in its "forking" path.
        """

        logger.debug("=====================")
        logger.debug("Setting main channels")
        logger.debug("=====================")

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

    def _set_secondary_inputs(self):
        """Sets the main raw inputs and secondary inputs on the init process

        This method will fetch the :class:`assemblerflow.process.Init` process
        instance and sets the raw input (
        :func:`assemblerflow.process.Init.set_raw_inputs`) and the secondary
        inputs (:func:`assemblerflow.process.Init.set_secondary_inputs`) for
        that process. This will handle the connection of the user parameters
        with channels that are then consumed in the pipeline.
        """

        logger.debug("========================")
        logger.debug("Setting secondary inputs")
        logger.debug("========================")

        # Get init process
        init_process = self.processes[0]
        logger.debug("Setting main raw inputs: "
                     "{}".format(self.main_raw_inputs))
        init_process.set_raw_inputs(self.main_raw_inputs)
        logger.debug("Setting secondary inputs: "
                     "{}".format(self.secondary_inputs))
        init_process.set_secondary_inputs(self.secondary_inputs)

    def _set_secondary_channels(self):
        """Sets the secondary channels for the pipeline

        This will iterate over the
        :py:attr:`NextflowGenerator.secondary_channels` dictionary that is
        populated when executing
        :func:`~NextflowGenerator._update_secondary_channels` method.
        """

        logger.debug("==========================")
        logger.debug("Setting secondary channels")
        logger.debug("==========================")

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

    def _set_compiler_channels(self):

        self._set_status_channels()

    def _set_status_channels(self):
        """Compiles all status channels for the status compiler process
        """

        status_inst = pc.StatusCompiler(template="status_compiler")

        # Compile status channels from pipeline process
        status_channels = []
        for p in [p for p in self.processes]:
            if not any([isinstance(p, x) for x in self.skip_class]):
                status_channels.extend(p.status_strs)

        if not status_channels:
            logger.debug("No status channels found. Skipping status compiler"
                         "process")
            return

        logger.debug("Setting status channels: {}".format(status_channels))

        # Check for duplicate channels. Raise exception if found.
        if len(status_channels) != len(set(status_channels)):
            raise eh.ProcessError(
                "Duplicate status channels detected. Please ensure that "
                "the 'status_channels' attributes of each process are "
                "unique. Here are the status channels:\n\n{}".format(
                    ", ".join(status_channels)
                ))

        status_inst.set_status_channels(status_channels)
        self.processes.append(status_inst)

    @staticmethod
    def _get_resources_string(res_dict, pid):
        """ Returns the nextflow resources string from a dictionary object

        If the dictionary has at least on of the resource directives, these
        will be compiled for each process in the dictionary and returned
        as a string read for injection in the nextflow config file template.

        This dictionary should be::

            dict = {"processA": {"cpus": 1, "memory": "4GB"},
                    "processB": {"cpus": 2}}

        Parameters
        ----------
        res_dict : dict
            Dictionary with the resources for processes.
        pid : int
            Unique identified of the process

        Returns
        -------
        str : nextflow config string
        """

        resource_directives = ["cpus", "memory"]
        config_str = ""

        for p, directives in res_dict.items():

            for d, val in directives.items():

                if d not in resource_directives:
                    continue

                if "{" in str(val):
                    config_str += '\n${}_{}.{} = {}'.format(p, pid, d, val)
                else:
                    config_str += '\n${}_{}.{} = "{}"'.format(p, pid, d, val)

        return config_str

    @staticmethod
    def _get_container_string(cont_dict, pid):
        """ Returns the nextflow containers string from a dictionary object

        If the dictionary has at least on of the container directives, these
        will be compiled for each process in the dictionary and returned
        as a string read for injection in the nextflow config file template.

        This dictionary should be::

            dict = {"processA": {"container": "asd", "version": "1.0.0"},
                    "processB": {"container": "dsd"}}

        Parameters
        ----------
        cont_dict : dict
            Dictionary with the containers for processes.
        pid : int
            Unique identified of the process

        Returns
        -------
        str : nextflow config string
        """

        config_str = ""

        for p, directives in cont_dict.items():

            container = ""

            if "container" in directives:
                container += directives["container"]

                if "version" in directives:
                    container += ":{}".format(directives["version"])
                else:
                    container += ":latest"

            config_str += '\n${}_{}.container = "{}"'.format(p, pid, container)

        return config_str

    @staticmethod
    def _render_config(template, context):

        tpl_dir = join(dirname(abspath(__file__)), "templates")
        tpl_path = join(tpl_dir, template + ".config")

        path, filename = split(tpl_path)

        return jinja2.Environment(
            loader=jinja2.FileSystemLoader(path or "./")
        ).get_template(filename).render(context)

    def _set_configurations(self):
        """This method will iterate over all process in the pipeline and
        populate the nextflow configuration files with the directives
        of each process in the pipeline.
        """

        resources = ""
        containers = ""

        for p in self.processes:

            # Skip processes with the directives attribute populated
            if not p.directives:
                continue

            resources += self._get_resources_string(p.directives, p.pid)
            containers += self._get_container_string(p.directives, p.pid)

        self.resources = self._render_config("resources", {
            "process_info": resources
        })
        self.containers = self._render_config("containers", {
            "container_info": containers
        })

    def build(self):
        """Main pipeline builder

        This method is responsible for building the
        :py:attr:`NextflowGenerator.template` attribute that will contain
        the nextflow code of the pipeline.

        First it builds the header, then sets the main channels, the
        secondary inputs, secondary channels and finally the
        status channels. When the pipeline is built, is writes the code
        to a nextflow file.
        """

        # Generate regular nextflow header that sets up the shebang, imports
        # and all possible initial channels
        self._build_header()

        self._set_channels()

        self._set_secondary_inputs()

        self._set_secondary_channels()

        self._set_compiler_channels()

        self._set_configurations()

        for p in self.processes:
            self.template += p.template_str

        self._build_footer()

        project_root = dirname(self.nf_file)

        # Write pipeline file
        with open(self.nf_file, "w") as fh:
            fh.write(self.template)

        # Write resources config
        with open(join(project_root, "resources.config"), "w") as fh:
            fh.write(self.resources)

        # Write containers config
        with open(join(project_root, "containers.config"), "w") as fh:
            fh.write(self.containers)
