import os
import jinja2
import logging

from os.path import dirname, join, abspath

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
    ptype : str
        Process type. See :py:attr:`Process.accepted_types`.
    template : str
        Name of the jinja2 template with the nextflow code for that process.
        Templates are stored in ``generator/templates``.
    """

    def __init__(self, ptype, template, process_id=None):

        accepted_types = [
            "init",
            "raw",
            "pre_assembly",
            "assembly",
            "post_assembly",
            "annotation",
            "status",
            "terminal"
        ]
        """
        list: Accepted process types
        """

        if ptype not in accepted_types:
            raise ValueError(
                "{} is not an accepted process type".format(ptype))

        self.pid = None
        """
        int: Process ID number that represents the order and position in the
        generated pipeline
        """

        self.process_id = process_id
        """
        int or str: optional Process ID that has no effect on the setup of
        the pipeline channels. It's used for the POST requests of each main
        process and is mapped to the process IDs of the innuendo/oneida
        platform
        """

        self.ptype = ptype
        """
        str: Process type. See :py:attr:`accepted_types`.
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

        self._main_in_str = None
        """
        str: String used to specify the prefix of main input channel.
        """

        self._main_out_str = None
        """
        str: String used to specify the prefix of main output channel.
        """

        self._input_channel = None
        """
        str: Place holder of the main input channel for the current process.
        This attribute can change dynamically depending on the forks and
        secondary channels in the final pipeline.
        """

        self._output_channel = None
        """
        str: Place holder of the main output channel for the current process.
        This attribute can change dynamically depending on the forks and
        secondary channels in the final pipeline.
        """

        self._set_main_channel_name(ptype)

        self.link_start = [self._main_out_str]
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

        self.status_channels = ["STATUS"]
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

        self._context = None
        """
        dict: Dictionary with the keyword placeholders for the string template
        of the current process.
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
            raise Exception("Template {} does not exist".format(tpl_path))

        self._template_path = join(tpl_dir, template + ".nf")

    def _set_main_channel_name(self, ptype):
        """Sets the prefix for the main channel depending on the process type

        ``Pre-assembly`` types are set to ``MAIN_fq``, while ``post-assembly``
        are set to ``MAIN_assembly``. This distinction is important to allow
        the forking of the last main channel with FastQ files or with
        assembly files.
        """

        if ptype == "init":
            self._main_in_str = "MAIN_raw"
            self._main_out_str = "MAIN_raw"
        elif ptype == "raw":
            self._main_in_str = self._main_out_str = "MAIN_raw"
        elif ptype == "pre_assembly":
            self._main_in_str = self._main_out_str = "MAIN_fq"
        elif ptype == "assembly":
            self._main_in_str = "MAIN_fq"
            self._main_out_str = "MAIN_assembly"
        else:
            self._main_in_str = self._main_out_str = "MAIN_assembly"

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
            raise Exception("Channels must be setup first using the "
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

        self.pid = kwargs.get("pid")

        self._input_channel = "{}_{}".format(self._main_in_str, self.pid)

        for i in self.status_channels:
            self.status_strs.append("{}_{}".format(i, self.pid))

        if self.output_type:
            self._output_channel = "{}_{}".format(self._main_out_str,
                                                  self.pid + 1)

        self._context = {**kwargs, **{"input_channel": self._input_channel,
                                      "output_channel": self._output_channel,
                                      "template": self.template}}

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

        # Handle the case where the main channel is forked
        if source.startswith("MAIN"):
            # Update previous output_channel to prevent overlap with
            # subsequent main channels. This is done by adding a "_" at the
            # beginning of the channel name
            self._context["output_channel"] = "_{}".format(
                self._output_channel)
            # Set source to modified output channel
            source = self._context["output_channel"]
            # Add the next first main channel to the channel_list
            channel_list.append(self._output_channel)
        # Handle forks from non main channels
        else:
            source = "{}_{}".format(source, self.pid)

        # Removes possible duplicate channels, when the fork is terminal
        channel_list = list(set(channel_list))

        # When there is only one channel to fork into, use the 'set' operator
        # instead of 'into'
        if len(channel_list) == 1:
            self.forks.append("\n{}.set{{ {} }}\n".format(source,
                                                           channel_list[0]))
        else:
            self.forks.append("\n{}.into{{ {} }}\n".format(
                source, ";".join(channel_list)))

        logger.debug("Setting forks attribute to: {}".format(self.forks))
        self._context = {**self._context, **{"forks": "\n".join(self.forks)}}


class Status(Process):
    """Extends the Process methods to status-type processes
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

    def set_status_channels(self, channel_list):
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
        """

        if len(channel_list) == 1:
            logger.debug("Setting only one status channel: {}".format(
                channel_list[0]))
            self._context = {"status_channels": channel_list[0]}

        else:
            first_status = channel_list[0]
            lst = ",".join(channel_list[1:])

            s = "{}.mix({})".format(first_status, lst)

            logger.debug("Status channel string: {}".format(s))

            self._context = {"status_channels": s}


class Init(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="init", **kwargs)

        self.input_type = None
        self.output_type = "raw"

        self.status_channels = []

    def set_secondary_channel(self, source, channel_list):

        logger.debug("Setting secondary channel for source '{}': {}".format(
            source, channel_list))

        if len(channel_list) == 1:
            self.forks.append("\nIN_fastq_raw.set{{ {} }}\n".format(
                channel_list[0]))
        else:
            self.forks.append("\nIN_fastq_raw.into{{ {} }}\n".format(
                ";".join(channel_list)
            ))

        logger.debug("Setting forks attribute to: {}".format(self.forks))
        self._context = {**self._context, **{"forks": "\n".join(self.forks)}}
        logger.debug(self._context)


class IntegrityCoverage(Process):
    """Process template interface for first integrity_coverage process

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains two **secondary channel link starts**:

        - ``SIDE_phred``: Phred score of the FastQ files
        - ``SIDE_max_len``: Maximum read length
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "raw"
        self.output_type = "fastq"

        self._main_in_str = "MAIN_raw"

        self.link_start.extend(["SIDE_phred", "SIDE_max_len"])

        self.link_end.append({"link": "MAIN_raw",
                              "alias": "MAIN_raw"})


class SeqTyping(Process):
    """

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="raw", **kwargs)

        self.input_type = "raw"
        self.output_type = None

        self.ignore_type = True
        self.ignore_pid = True

        self.status_channels = []

        self.link_start = None
        self.link_end.append({"link": "MAIN_raw",
                              "alias": "SIDE_SeqType_raw"})


class PathoTyping(Process):
    """

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="raw", **kwargs)

        self.input_type = "raw"
        self.output_type = None

        self.ignore_type = True
        self.ignore_pid = True

        self.status_channels = []

        self.link_start = None
        self.link_end.append({"link": "MAIN_raw",
                              "alias": "SIDE_PathoType_raw"})


class CheckCoverage(Process):
    """Process template interface for additional integrity_coverage process

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains one **secondary channel link start**:

        - ``SIDE_max_len``: Maximum read length

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_start.extend(["SIDE_max_len"])


class FastQC(Process):
    """FastQC process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains two **status channels**:

        - ``STATUS_fastqc``: Status for the fastqc process
        - ``STATUS_report``: Status for the fastqc_report process

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.status_channels = ["STATUS_fastqc", "STATUS_report"]
        """
        list: Setting status channels for FastQC execution and FastQC report
        """


class Trimmomatic(Process):
    """Trimmomatic process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains one **secondary channel link end**:

        - ``SIDE_phred`` (alias: ``SIDE_phred``): Receives FastQ phred score
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_end.append({"link": "SIDE_phred", "alias": "SIDE_phred"})


class FastqcTrimmomatic(Process):
    """Fastqc + Trimmomatic process template interface

    This process executes FastQC only to inform the trim range for trimmomatic,
    not for QC checks.

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains one **secondary channel link end**:

        - ``SIDE_phred`` (alias: ``SIDE_phred``): Receives FastQ phred score

    It contains three **status channels**:

        - ``STATUS_fastqc``: Status for the fastqc process
        - ``STATUS_report``: Status for the fastqc_report process
        - ``STATUS_trim``: Status for the trimmomatic process
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_end.append({"link": "SIDE_phred", "alias": "SIDE_phred"})

        self.status_channels = ["STATUS_fastqc", "STATUS_report",
                                "STATUS_trim"]


class Spades(Process):
    """Spades process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: assembly
        - ``ptype``: assembly

    It contains one **secondary channel link end**:

        - ``SIDE_max_len`` (alias: ``SIDE_max_len``): Receives max read length
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "assembly"

        self.link_end.append({"link": "SIDE_max_len", "alias": "SIDE_max_len"})


class ProcessSpades(Process):
    """Process spades process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: assembly
        - ``ptype``: post_assembly

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="post_assembly",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = "assembly"


class AssemblyMapping(Process):
    """Assembly mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: assembly
        - ``ptype``: post_assembly

    It contains one **secondary channel link end**:

        - ``MAIN_fq`` (alias: ``_MAIN_assembly``): Receives the FastQ files
        from the last process with ``fastq`` output type.

    It contains two **status channels**:

        - ``STATUS_am``: Status for the assembly_mapping process
        - ``STATUS_amp``: Status for the process_assembly_mapping process
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="post_assembly",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = "assembly"

        self.status_channels = ["STATUS_am", "STATUS_amp"]

        self.link_start.append("SIDE_BpCoverage")
        self.link_end.append({"link": "MAIN_fq", "alias": "_MAIN_assembly"})


class Pilon(Process):
    """Pilon mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: assembly
        - ``ptype``: post_assembly

    It contains one **dependency process**:

        - ``assembly_mapping``: Requires the BAM file generated by the
        assembly mapping process
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="post_assembly",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = "assembly"

        self.dependencies = ["assembly_mapping"]

        self.link_end.append({"link": "SIDE_BpCoverage",
                              "alias": "SIDE_BpCoverage"})


class Mlst(Process):
    """Mlst mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: None
        - ``ptype``: post_assembly

    It contains one **secondary channel link end**:

        - ``MAIN_assembly`` (alias: ``MAIN_assembly``): Receives the last
        assembly.
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="annotation",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = "assembly"


class Abricate(Process):
    """Abricate mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: None
        - ``ptype``: post_assembly

    It contains one **secondary channel link end**:

        - ``MAIN_assembly`` (alias: ``MAIN_assembly``): Receives the last
        assembly.
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="annotation",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = None

        self.ignore_type = True

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})


class Prokka(Process):
    """Prokka mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: None
        - ``ptype``: post_assembly

    It contains one **secondary channel link end**:

        - ``MAIN_assembly`` (alias: ``MAIN_assembly``): Receives the last
        assembly.
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="annotation",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = None

        self.ignore_type = True

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})


class Chewbbaca(Process):
    """Prokka mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: None
        - ``ptype``: post_assembly

    It contains one **secondary channel link end**:

        - ``MAIN_assembly`` (alias: ``MAIN_assembly``): Receives the last
        assembly.
    """

    def __init__(self, **kwargs):

        super().__init__(ptype="annotation",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = None

        self.ignore_type = True

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})


class TraceCompiler(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="terminal", **kwargs)

        self.link_start = None

        self.ignore_type = True


class StatusCompiler(Status):
    """Status compiler process template interface

    This special process receives the status channels from all processes
    in the generated pipeline.

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="status",
                         **kwargs)

        self.ignore_type = True

        self.link_start = None

