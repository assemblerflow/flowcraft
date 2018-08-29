import os
import jinja2
import logging

from os.path import dirname, join, abspath

try:
    import generator.error_handling as eh
except ImportError:
    import flowcraft.generator.error_handling as eh

logger = logging.getLogger("main.{}".format(__name__))


class Process:
    """Main interface for basic process functionality

    The ``Process`` class is intended to be inherited by specific process
    classes (e.g., :py:class:`IntegrityCoverage`) and provides the basic
    functionality to build the channels and links between processes.

    Child classes are expected to inherit the ``__init__`` execution, which
    basically means that at least, the child must be defined as::

        class ChildProcess(Process):
            def__init__(self, **kwargs):
                super().__init__(**kwargs)

    This ensures that when the ``ChildProcess`` class is instantiated, it
    automatically sets the attributes of the parent class.

    This also means that child processes must be instantiated providing
    information on the process type and jinja2 template with the nextflow code.

    Parameters
    ----------
    template : str
        Name of the jinja2 template with the nextflow code for that process.
        Templates are stored in ``generator/templates``.
    """

    RAW_MAPPING = {
        "fastq": {
            "params": "fastq",
            "description": "Path expression to paired-end fastq files."
                           " (default: $params.fastq)",
            "default_value": "'fastq/*_{1,2}.*'",
            "channel": "IN_fastq_raw",
            "channel_str":
                "Channel.fromFilePairs(params.{0})"
                ".ifEmpty {{ exit 1, \"No fastq files provided with pattern:"
                "'${{params.{0}}}'\" }}",
            "checks":
                "if (params.{0} instanceof Boolean){{"
                "exit 1, \"'{0}' must be a path pattern. Provide value:"
                "'$params.{0}'\"}}\n"
                "if (!params.{0}){{ exit 1, \"'{0}' parameter "
                "missing\"}}"
        },
        "fasta": {
            "params": "fasta",
            "description": "Path fasta files. (default: $params.fastq)",
            "default_value": "'fasta/*.fasta'",
            "channel": "IN_fasta_raw",
            "channel_str":
                "Channel.fromPath(params.{0})."
                "map{{ it -> file(it).exists() ? [it.toString()"
                ".tokenize('/').last()"
                ".tokenize('.')[0..-2].join('.'), it] : null }}"
                ".ifEmpty {{ exit 1, \"No fasta files provided with pattern:"
                "'${{params.{0}}}'\" }}",
            "checks":
                "if (params.{0} instanceof Boolean){{"
                "exit 1, \"'{0}' must be a path pattern. Provide value:"
                "'$params.{0}'\"}}\n"
                "if (!params.{0}){{ exit 1, \"'{0}' parameter "
                "missing\"}}"
        },
        "accessions": {
            "params": "accessions",
            "description": "Path file with accessions, one perline. ("
                           "default: $params.fastq)",
            "default_value": "null",
            "channel": "IN_accessions_raw",
            "channel_str":
                "Channel.fromPath(params.{0})"
                ".ifEmpty {{ exit 1, \"No accessions file provided with path:"
                "'${{params.{0}}}'\" }}",
            "checks":
                "if (!params.{0}){{ exit 1, \"'{0}' parameter "
                "missing\" }}\n"
        }
    }
    """
    dict: Contains the mapping between the :attr:`Process.input_type` attribute
    and the corresponding nextflow parameter and main channel definition,
    e.g.::

        "fastq" : {
            "params": "fastq",
            "channel: "<channel>
        }
    """

    def __init__(self, template):

        self.pid = None
        """
        int: Process ID number that represents the order and position in the
        generated pipeline
        """

        self.template = template
        """
        str: Template name for the current process. This string will be used
        to fetch the file containing the corresponding jinja2 template
        in the :py:func:`_set_template` method
        """

        self._template_path = None
        """
        str: Path to the file containing the jinja2 template file. It's
        set in :py:func:`_set_template`.
        """
        self._set_template(template)

        self.input_type = None
        """
        str: Type of expected input data. Used to verify the connection between
        two processes is viable.
        """

        self.output_type = None
        """
        str: Type of output data. Used to verify the connection between
        two processes is viable.
        """

        self.ignore_type = False
        """
        boolean: If True, this process will ignore the input/output type
        requirements. This attribute is set to True for terminal singleton
        forks in the pipeline.
        """

        self.ignore_pid = False
        """
        boolean: If True, this process will not make the pid advance. This
        is used for terminal forks before the end of the pipeline.
        """

        self.dependencies = []
        """
        list: Contains the dependencies of the current process in the form
        of the :py:attr:`Process.template` attribute (e.g., [``fastqc``])
        """

        self.lane = None
        self.parent_lane = None

        self.input_channel = None
        """
        str: Place holder of the main input channel for the current process.
        This attribute can change dynamically depending on the forks and
        secondary channels in the final pipeline.
        """

        self.output_channel = None
        """
        str: Place holder of the main output channel for the current process.
        This attribute can change dynamically depending on the forks and
        secondary channels in the final pipeline.
        """

        self.input_user_channel = None
        """
        dict: Stores a dictionary of two key:value pairs containing
        the raw input channel for the process. This is automatically
         determined by the :attr:`~Process.input_type` attribute, and will
        fetch the information that is mapped in the :attr:`RAW_MAPPING`
         variable. It will only be used by the first process(es) defined in
         a pipeline.
        """

        self.link_start = []
        """
        list: List of strings with the starting points for secondary channels.
        When building the pipeline, these strings will be matched with equal
        strings in the :py:attr:`link_end` attribute of other Processes.
        """

        self.link_end = []
        """
        list: List of dictionaries containing the a string of the ending point
        for a secondary channel. Each dictionary should contain at least
        two key/vals:
        ``{"link": <link string>, "alias":<string for template>}``
        """

        self.status_channels = ["STATUS_{}".format(template)]
        """
        list: Name of the status channels produced by the process. By default,
        it sets a single status channel. If more than one status channels
        are required for the process, list each one in this attribute
        (e.g., :py:attr:`FastQC.status_channels`)
        """
        self.status_strs = []
        """
        str: Name of the status channel for the current process. These strings
        will be provided to the StatusCompiler process to collect and
        compile status reports
        """

        self.forks = []
        """
        list: List of strings with the literal definition of the forks for
        the current process, ready to be added to the template string.
        """

        self.main_forks = []
        """
        list: List of the channels onto which the main output should be
        forked into. They will be automatically added to the
        :attr:`~Process.main_forks` attribute when setting the secondary
        channels
        """

        self.secondary_inputs = []
        """
        list: List of dictionaries with secondary input channels from nextflow
        parameters. This dictionary should contain two key:value pairs
        with the ``params`` key, containing the parameter name, and the
        ``channel`` key, containing the nextflow channel definition::

            {
                "params": "pathoSpecies",
                "channel": "IN_pathoSpecies = Channel
                                                .value(params.pathoSpecies)"
            }
        """
        self.secondary_input_str = ""

        self.extra_input = ""
        """
        str:  with the name of the params that will be used to provide
        extra input into the process. This extra input will be mixed with
        the main input channel using nextflow's ``mix`` operator. Its
        channel will be defined at the start of the pipeline, based on the
        ``channel_str`` key of the :attr:`~Process.RAW_MAPPING` for the
        corresponding input type.
        """

        self.params = {}
        """
        dict: Maps the parameter names to the corresponding default values.
        """

        self.param_id = ""
        """
        str: The parameter id suffix that will be added to each parameter. In
        case it is empty, the multiple identical parameters in different
        components will be merged.
        """

        self._context = {}
        """
        dict: Dictionary with the keyword placeholders for the string template
        of the current process.
        """

        self.directives = {}
        """
        dict: Specifies the directives (cpus, memory, container) for each
        nextflow process in the template. If specified, this directives
        will be added to the nextflow configuration file. Otherwise,
        the default values for cpus and memory will be used. In the case
        of containers, they will not run inside any container.

        The current supported directives are:
            - cpus
            - memory
            - container
            - container tag/version

        An example of directives for two process is as follows::
        
            self.directives = {
                "processA": {"cpus": 1, "memory": "1GB"},
                "processB": {"memory": "5GB", "container": "my/image",
                             "version": "0.5.0"}
            }
        """

        self.compiler = {}
        """
        dict: Specifies channels from the current process that are received
        by a compiler process. Each key in this dictionary should match
        a compiler process key in
        :attr:`~flowcraft.generator.engine.NextflowGenerator.compilers`.
        The value should be a list of the channels that will be fed to the
        compiler process::
        
            self.compiler["patlas_consensus"] = ["mashScreenOutputChannel"]
        """

    def _set_template(self, template):
        """Sets the path to the appropriate jinja template file

        When a Process instance is initialized, this method will fetch
        the location of the appropriate template file, based on the
        ``template`` argument. It will raise an exception is the template
        file is not found. Otherwise, it will set the
        :py:attr:`Process.template_path` attribute.
        """

        # Set template directory
        tpl_dir = join(dirname(abspath(__file__)), "templates")

        # Set template file path
        tpl_path = join(tpl_dir, template + ".nf")

        if not os.path.exists(tpl_path):
            raise eh.ProcessError(
                "Template {} does not exist".format(tpl_path))

        self._template_path = join(tpl_dir, template + ".nf")

    def set_main_channel_names(self, input_suffix, output_suffix, lane):
        """Sets the main channel names based on the provide input and
        output channel suffixes. This is performed when connecting processes.

        Parameters
        ----------
        input_suffix : str
            Suffix added to the input channel. Should be based on the lane
            and an arbitrary unique id
        output_suffix : str
            Suffix added to the output channel. Should be based on the lane
            and an arbitrary unique id
        lane : int
            Sets the lane of the process.
        """

        self.input_channel = "{}_in_{}".format(self.template, input_suffix)
        self.output_channel = "{}_out_{}".format(self.template, output_suffix)
        self.lane = lane

    def set_param_id(self, param_id):
        """Sets the param_id for the process, which will be used to render
        the template.

        Parameters
        ----------
        param_id : str
            The :attr:`param_id` attribute of the process.
        """

        self._context = {**self._context, "param_id": param_id}

    def get_user_channel(self, input_channel, input_type=None):
        """Returns the main raw channel for the process

        Provided with at least a channel name, this method returns the raw
        channel name and specification (the nextflow string definition)
        for the process. By default, it will fork from the raw input of
        the process' :attr:`~Process.input_type` attribute. However, this
        behaviour can be overridden by providing the ``input_type`` argument.

        If the specified or inferred input type exists in the
        :attr:`~Process.RAW_MAPPING` dictionary, the channel info dictionary
        will be retrieved along with the specified input channel. Otherwise,
        it will return None.

        An example of the returned dictionary is::

             {"input_channel": "myChannel",
             "params": "fastq",
             "channel": "IN_fastq_raw",
             "channel_str":"IN_fastq_raw = Channel.fromFilePairs(params.fastq)"
            }

        Returns
        -------
        dict or None
            Dictionary with the complete raw channel info. None if no
            channel is found.
        """

        res = {"input_channel": input_channel}

        itype = input_type if input_type else self.input_type

        if itype in self.RAW_MAPPING:

            channel_info = self.RAW_MAPPING[itype]

            return {**res, **channel_info}

    @staticmethod
    def render(template, context):
        """Wrapper to the jinja2 render method from a template file

        Parameters
        ----------
        template : str
            Path to template file.
        context : dict
            Dictionary with kwargs context to populate the template
        """

        path, filename = os.path.split(template)

        return jinja2.Environment(
            loader=jinja2.FileSystemLoader(path or './')
        ).get_template(filename).render(context)

    @property
    def template_str(self):
        """Class property that returns a populated template string

        This property allows the template of a particular process to be
        dynamically generated and returned when doing ``Process.template_str``.

        Returns
        -------
        x : str
            String with the complete and populated process template

        """

        if not self._context:
            raise eh.ProcessError("Channels must be setup first using the "
                                  "set_channels method")

        logger.debug("Setting context for template {}: {}".format(
            self.template, self._context
        ))

        x = self.render(self._template_path, self._context)
        return x

    def set_channels(self, **kwargs):
        """ General purpose method that sets the main channels

        This method will take a variable number of keyword arguments to
        set the :py:attr:`Process._context` attribute with the information
        on the main channels for the process. This is done by appending
        the process ID (:py:attr:`Process.pid`) attribute to the input,
        output and status channel prefix strings. In the output channel,
        the process ID is incremented by 1 to allow the connection with the
        channel in the next process.

        The ``**kwargs`` system for setting the :py:attr:`Process._context`
        attribute also provides additional flexibility. In this way,
        individual processes can provide additional information not covered
        in this method, without changing it.

        Parameters
        ----------
        kwargs : dict
            Dictionary with the keyword arguments for setting up the template
            context
        """

        if not self.pid:
            self.pid = "{}_{}".format(self.lane, kwargs.get("pid"))

        for i in self.status_channels:
            if i.startswith("STATUS_"):
                self.status_strs.append("{}_{}".format(i, self.pid))
            else:
                self.status_strs.append("STATUS_{}_{}".format(i, self.pid))

        if self.main_forks:
            logger.debug("Setting main fork channels: {}".format(
                self.main_forks))
            operator = "set" if len(self.main_forks) == 1 else "into"
            self.forks = ["\n{}.{}{{ {} }}\n".format(
                self.output_channel, operator, ";".join(self.main_forks))]

        self._context = {**kwargs, **{"input_channel": self.input_channel,
                                      "output_channel": self.output_channel,
                                      "template": self.template,
                                      "forks": "\n".join(self.forks),
                                      "pid": self.pid}}

    def update_main_input(self, input_str):

        self.input_channel = input_str
        self._context["input_channel"] = self.input_channel

    def update_main_forks(self, sink):
        """Updates the forks attribute with the sink channel destination

        Parameters
        ----------
        sink : str
            Channel onto which the main input will be forked to

        """

        if not self.main_forks:
            self.main_forks = [self.output_channel]
            self.output_channel = "_{}".format(self.output_channel)
        self.main_forks.append(sink)

        # fork_lst = self.forks + self.main_forks
        operator = "set" if len(self.main_forks) == 1 else "into"
        self.forks = ["\n{}.{}{{ {} }}\n".format(
            self.output_channel, operator, ";".join(self.main_forks))]

        self._context = {**self._context,
                         **{"forks": "".join(self.forks),
                            "output_channel": self.output_channel}}

    def set_secondary_channel(self, source, channel_list):
        """ General purpose method for setting a secondary channel

        This method allows a given source channel to be forked into one or
        more channels and sets those forks in the :py:attr:`Process.forks`
        attribute. Both the source and the channels in the ``channel_list``
        argument must be the final channel strings,  which means that this
        method should be called only after setting the main channels.

        If the source is not a main channel, this will simply create a fork
        or set for every channel in the ``channel_list`` argument list::

            SOURCE_CHANNEL_1.into{SINK_1;SINK_2}

        If the source is a main channel, this will apply some changes to
        the output channel of the process, to avoid overlapping main output
        channels.  For instance, forking the main output channel for process
        2 would create a ``MAIN_2.into{...}``. The issue here is that the
        ``MAIN_2`` channel is expected as the input of the next process, but
        now is being used to create the fork. To solve this issue, the output
        channel is modified into ``_MAIN_2``, and the fork is set to
        the channels provided channels plus the ``MAIN_2`` channel::

            _MAIN_2.into{MAIN_2;MAIN_5;...}

        Parameters
        ----------
        source : str
            String with the name of the source channel
        channel_list : list
            List of channels that will receive a fork of the secondary
            channel
        """

        logger.debug("Setting secondary channel for source '{}': {}".format(
            source, channel_list))

        source = "{}_{}".format(source, self.pid)

        # Removes possible duplicate channels, when the fork is terminal
        channel_list = sorted(list(set(channel_list)))

        # When there is only one channel to fork into, use the 'set' operator
        # instead of 'into'
        op = "set" if len(channel_list) == 1 else "into"
        self.forks.append("\n{}.{}{{ {} }}\n".format(
            source, op, ";".join(channel_list)))

        logger.debug("Setting forks attribute to: {}".format(self.forks))
        self._context = {**self._context, **{"forks": "\n".join(self.forks)}}

    def update_attributes(self, attr_dict):
        """Updates the directives attribute from a dictionary object.

        This will only update the directives for processes that have been
        defined in the subclass.

        Parameters
        ----------
        attr_dict : dict
            Dictionary containing the attributes that will be used to update
            the process attributes and/or directives.

        """

        # Update directives
        # Allowed attributes to write
        valid_directives = ["pid", "ignore_type", "ignore_pid", "extra_input",
                            "group"]

        for attribute, val in attr_dict.items():

            # If the attribute has a valid directive key, update that
            # directive
            if attribute in valid_directives and hasattr(self, attribute):
                setattr(self, attribute, val)

            else:
                for p in self.directives:
                    self.directives[p][attribute] = val


class Compiler(Process):
    """Extends the Process methods to status-type processes
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.ignore_type = True
        self.link_start = None

    def set_compiler_channels(self, channel_list, operator="mix"):
        """General method for setting the input channels for the status process

        Given a list of status channels that are gathered during the pipeline
        construction, this method will automatically set the input channel
        for the status process. This makes use of the ``mix`` channel operator
        of nextflow for multiple channels::

            STATUS_1.mix(STATUS_2,STATUS_3,...)

        This will set the ``status_channels`` key for the ``_context``
        attribute of the process.

        Parameters
        ----------
        channel_list : list
            List of strings with the final name of the status channels
        operator : str
            Specifies the operator used to join the compiler channels.
            Available options are 'mix'and 'join'.
        """

        if not channel_list:
            raise eh.ProcessError("At least one status channel must be "
                                  "provided to include this process in the "
                                  "pipeline")

        if len(channel_list) == 1:
            logger.debug("Setting only one status channel: {}".format(
                channel_list[0]))
            self._context = {"compile_channels": channel_list[0]}

        else:

            first_status = channel_list[0]

            if operator == "mix":
                lst = ",".join(channel_list[1:])

                s = "{}.mix({})".format(first_status, lst)

            elif operator == "join":

                s = first_status
                for ch in channel_list[1:]:
                    s += ".join({})".format(ch)

                s += ".map{ ot -> [ ot[0], ot[1..-1] ] }"

            logger.debug("Status channel string: {}".format(s))

            self._context = {"compile_channels": s}


class Init(Process):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = None
        self.output_type = "raw"

        self.status_channels = []

    def set_raw_inputs(self, raw_input):
        """Sets the main input channels of the pipeline and their forks.

        The ``raw_input`` dictionary input should contain one entry for each
        input type (fastq, fasta, etc). The corresponding value should be a
        dictionary/json with the following key:values:

        - ``channel``: Name of the raw input channel (e.g.: channel1)
        - ``channel_str``: The nextflow definition of the channel and
           eventual checks (e.g.: channel1 = Channel.fromPath(param))
        - ``raw_forks``: A list of channels to which the channel name will
          for to.

        Each new type of input parameter is automatically added to the
        :attr:`params` attribute, so that they are automatically collected
        for the pipeline description and help.

        Parameters
        ----------
        raw_input : dict
            Contains an entry for each input type with the channel name,
            channel string and forks.
        """

        logger.debug("Setting raw inputs using raw input dict: {}".format(
            raw_input))

        primary_inputs = []

        for input_type, el in raw_input.items():

            primary_inputs.append(el["channel_str"])

            # Update the process' parameters with the raw input
            raw_channel = self.RAW_MAPPING[input_type]
            self.params[input_type] = {
                "default": raw_channel["default_value"],
                "description": raw_channel["description"]
            }

            op = "set" if len(el["raw_forks"]) == 1 else "into"

            self.forks.append("\n{}.{}{{ {} }}\n".format(
                el["channel"], op, ";".join(el["raw_forks"])
            ))

        logger.debug("Setting raw inputs: {}".format(primary_inputs))
        logger.debug("Setting forks attribute to: {}".format(self.forks))
        self._context = {**self._context,
                         **{"forks": "\n".join(self.forks),
                            "main_inputs": "\n".join(primary_inputs)}}

    def set_secondary_inputs(self, channel_dict):
        """ Adds secondary inputs to the start of the pipeline.

        This channels are inserted into the pipeline file as they are
        provided in the values of the argument.

        Parameters
        ----------
        channel_dict : dict
            Each entry should be <parameter>: <channel string>.
        """

        logger.debug("Setting secondary inputs: {}".format(channel_dict))

        secondary_input_str = "\n".join(list(channel_dict.values()))
        self._context = {**self._context,
                         **{"secondary_inputs": secondary_input_str}}

    def set_extra_inputs(self, channel_dict):
        """Sets the initial definition of the extra input channels.

        The ``channel_dict`` argument should contain the input type and
        destination channel of each parameter (which is the key)::

            channel_dict = {
                "param1": {
                    "input_type": "fasta"
                    "channels": ["abricate_2_3", "chewbbaca_3_4"]
                }
            }

        Parameters
        ----------
        channel_dict : dict
            Dictionary with the extra_input parameter as key, and a dictionary
            as a value with the input_type and destination channels
        """

        extra_inputs = []

        for param, info in channel_dict.items():

            # Update the process' parameters with the raw input
            raw_channel = self.RAW_MAPPING[info["input_type"]]
            self.params[param] = {
                "default": raw_channel["default_value"],
                "description": raw_channel["description"]
            }

            channel_name = "IN_{}_extraInput".format(param)
            channel_str = self.RAW_MAPPING[info["input_type"]]["channel_str"]
            extra_inputs.append("{} = {}".format(channel_name,
                                                 channel_str.format(param)))

            op = "set" if len(info["channels"]) == 1 else "into"
            extra_inputs.append("{}.{}{{ {} }}".format(
                channel_name, op, ";".join(info["channels"])))

        self._context = {
            **self._context,
            **{"extra_inputs": "\n".join(extra_inputs)}
        }


class StatusCompiler(Compiler):
    """Status compiler process template interface

    This special process receives the status channels from all processes
    in the generated pipeline.
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)


class ReportCompiler(Compiler):
    """Reports compiler process template interface

    This special process receives the report channels from all processes
    in the generated pipeline.
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)


class PatlasConsensus(Compiler):
    """Patlas consensus compiler process template interface

    This special process receives the channels associated with the
    ``patlas_consensus`` key.
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)
