import re
import os
import sys
import curses
import signal
import locale
import socket
import logging
import hashlib
import requests
import json

from pympler import asizeof
from os.path import join, abspath
from time import gmtime, strftime, sleep
from collections import defaultdict, OrderedDict

try:
    import generator.error_handling as eh
    from generator.process_details import colored_print
except ImportError:
    import flowcraft.generator.error_handling as eh
    from flowcraft.generator.process_details import colored_print

locale.setlocale(locale.LC_ALL, "")
code = locale.getpreferredencoding()

logger = logging.getLogger("main.{}".format(__name__))


def signal_handler(screen):
    """This function is bound to the SIGINT signal (like ctrl+c) to graciously
    exit the program and reset the curses options.
    """

    if screen:
        screen.clear()
        screen.refresh()

        curses.nocbreak()
        screen.keypad(0)
        curses.echo()
        curses.endwin()

    print("Exiting flowcraft inspection... Bye")
    sys.exit(0)


class NextflowInspector:

    MAX_RETRIES = 1000

    def __init__(self, trace_file, refresh_rate, pretty=False, ip_addr=None):

        self.trace_file = trace_file
        """
        str: Path to nextflow trace file.
        """

        self.trace_sizestamp = None
        """
        str: Stores the sizestamp of the last modification of the trace file.
        This is used to parse the file only when it has changed.
        """

        self.refresh_rate = refresh_rate
        """
        float: Frequency (in seconds) that the curses screen will be refreshed.
        """

        self.stored_ids = []
        """
        list: Stores the task_ids that have already been parsed. It is used
        to skip them when parsing the trace files multiple times.
        """

        self.stored_log_ids = []
        """
        list: Stores the time stamps of the log file lines that were already
        parsed. It is used to skip parsing the log files multilpe times
        """

        self.trace_info = defaultdict(list)
        """
        dict: Main object that stores the status information for each process
        name in the trace file.
        """

        self.process_stats = {}
        """
        dict: Contains some statistics for each process.
        """

        self.processes = OrderedDict()
        """
        dict: Dictionary of processes from the pipeline with the status of the
        channel as the value. This information is retrieved from the
        .nextflow.log file in the :func:`_parser_pipeline_processes` method
        and updated in the :func:`_update_barrier_status` and
        :func:`_update_process_stats` and :func:`_update_submission_status`.
        """

        self.process_tags = {}
        """
        dict: Dictionary of processes with summary information for each tag
        it processes
        """

        self.samples = []
        """
        list: List of samples inferred from the pipeline.
        """

        self.skip_processes = ["status", "compile_status", "report",
                               "compile_reports", "fullConsensus",
                               "compile_status_buffer"]
        """
        list: List of special processes that should be skipped for inspection
        purposes.
        """

        self.log_file = ".nextflow.log"
        """
        str: Name of the nextflow log file.
        """

        self.log_sizestamp = None
        """
        str: Stores the sizestamp of the last modification of the nextflow
        log file. This is used to parse the file only when it has changed.
        """

        self.pipeline_tag = ""
        """
        str: Tag of the pipeline, parsed from .nextflow.log
        """

        self.log_retry = 0
        """
        int: Each time the log file is not found, this counter is
        increased. Only when it matches the :attr:`MAX_RETRIES` attribute
        does it raises a FileNotFoundError.
        """

        self.trace_retry = 0
        """
        int: Each time the log file is not found, this counter is 
        increased. Only when it matches the :attr:`MAX_RETRIES` attribute
        does it raises a FileNotFoundError.
        """

        self.pipeline_name = ""
        """
        str: Name of the nextflow pipeline file.
        """

        self.time_start = None
        """
        datetime.time object with the starting time of the pipeline.
        """

        self.time_stop = None
        """
        datetime.time object with the finish time of the pipeline. This
        attribute is only set when the pipeline is not running.
        """

        self.workdir = os.getcwd()
        """
        str: Path to the pipeline work directory
        """

        self.execution_command = None
        """
        str: The command used to execute the pipeline
        """

        self.nextflow_version = None
        """
        str: Nextflow's version string, as retrieved from the log file.
        """

        self.run_status = ""
        """
        str: Status of the pipeline. Can be either 'running', 'aborted',
        'error', 'complete'.
        """

        self.abort_cause = None
        """
        str or None: When :attr:`run_status` is "aborted", this attribute
        will contain the reason provided in the nextflow log. When this
        attribute is not None, it will also trigger the sending of the
        final lines of the nextflow log to broadcast.
        """

        if not ip_addr:
            self.app_address = "http://192.92.149.169:80/"
        else:
            self.app_address = ip_addr
            """
            str: Address of flowcraft web app
            """

        self.broadcast_address = "{}inspect/api/status".format(
            self.app_address)
        """
        str: Address of the REST api where the information will be sent
        """

        self._c = 0
        """
        Counter of payloads sent, for debug purposes
        """

        self.send = True
        """
        boolean: This attribute will be set to False after sending a request
        and set to True when there is a change in the inspection attributes.
        """

        # Skip these process names (they are check with the startswith()
        # method) when using the --pretty option
        if pretty:
            self._blacklist = [
                "report_coverage_", "fastqc2_report", "compile_fastqc_status2",
                "fastqc_report", "trim_report", "compile_fastqc_status",
                "report_corrupt_", "jsonDumpingMapping", "compile_mlst_",
                "mashOutputJson_", "mashDistOutputJson_", "pilon_report_",
                "compile_pilon_report"
            ]
        else:
            self._blacklist = []

        # CURSES ATTRIBUTES
        # Init curses screen
        self.screen = None
        self.top_line = 0
        self.padding = 0
        self.screen_lines = None
        self.max_width = 0
        self.content_lines = 0

        # Checks if nextflow log and trace files are available
        self._check_required_files()
        # Gathers the complete list of processes from the nextflow log
        self._get_pipeline_processes()
        # Fetches the pipeline status from the nextflow log
        self._update_pipeline_status()

        # Bind SIGINT to singal_handler function. This makes a clean exit
        # from the curses interface when exiting through ctrl+c.
        signal.signal(signal.SIGINT, lambda *x: signal_handler(self.screen))

    #################
    # UTILITY METHODS
    #################

    def _check_required_files(self):
        """Checks whetner the trace and log files are available
        """

        if not os.path.exists(self.trace_file):
            raise eh.InspectionError("The provided trace file could not be "
                                     "opened: {}".format(self.trace_file))

        if not os.path.exists(self.log_file):
            raise eh.InspectionError("The .nextflow.log files could not be "
                                     "opened. Are you sure you are in a "
                                     "nextflow project directory?")

    @staticmethod
    def _header_mapping(header):
        """Parses the trace file header and retrieves the positions of each
        column key.

        Parameters
        ----------
        header : str
            The header line of nextflow's trace file

        Returns
        -------
        dict
            Mapping the column ID to its position (e.g.: {"tag":2})
        """

        return dict(
            (x.strip(), pos) for pos, x in enumerate(header.split("\t"))
        )

    @staticmethod
    def _expand_path(hash_str):
        """Expands the hash string of a process (ae/1dasjdm) into a full
        working directory

        Parameters
        ----------
        hash_str : str
            Nextflow process hash with the beggining of the work directory

        Returns
        -------
        str
            Path to working directory of the hash string
        """

        try:
            first_hash, second_hash = hash_str.split("/")
            first_hash_path = join(abspath("work"), first_hash)

            for l in os.listdir(first_hash_path):
                if l.startswith(second_hash):
                    return join(first_hash_path, l)
        except FileNotFoundError:
            return None

    @staticmethod
    def _hms(s):
        """Converts a hms string into seconds.

        Parameters
        ----------
        s : str
            The hms string can be something like '20s', '1m30s' or '300ms'.

        Returns
        -------
        float
            Time in seconds.

        """

        if s == "-":
            return 0

        if s.endswith("ms"):
            return float(s.rstrip("ms")) / 1000

        fields = list(map(float, re.split("[hms]", s)[:-1]))
        if len(fields) == 3:
            return fields[0] * 3600 + fields[1] * 60 + fields[2]
        elif len(fields) == 2:
            return fields[0] * 60 + fields[1]
        else:
            return fields[0]

    @staticmethod
    def _size_coverter(s):
        """Converts size string into megabytes

        Parameters
        ----------
        s : str
            The size string can be '30KB', '20MB' or '1GB'

        Returns
        -------
        float
            With the size in bytes

        """

        if s.upper().endswith("KB"):
            return float(s.rstrip("KB")) / 1024

        elif s.upper().endswith(" B"):
            return float(s.rstrip("B")) / 1024 / 1024

        elif s.upper().endswith("MB"):
            return float(s.rstrip("MB"))

        elif s.upper().endswith("GB"):
            return float(s.rstrip("GB")) * 1024

        elif s.upper().endswith("TB"):
            return float(s.rstrip("TB")) * 1024 * 1024

        else:
            return float(s)

    @staticmethod
    def _size_compress(s):
        """Shortens a megabytes string.
        """

        if s / 1024 > 1:
            return "{}GB".format(round(s / 1024, 1))
        else:
            return "{}MB".format(s)

    #########################
    # AUXILIARY PARSE METHODS
    #########################

    def _get_pipeline_processes(self):
        """Parses the .nextflow.log file and retrieves the complete list
        of processes

        This method searches for specific signatures at the beginning of the
        .nextflow.log file::

             Apr-19 19:07:32.660 [main] DEBUG nextflow.processor
             TaskProcessor - Creating operator > report_corrupt_1_1 --
             maxForks: 4

        When a line with the .*Creating operator.* signature is found, the
        process name is retrieved and populates the :attr:`processes` attribute
        """

        with open(self.log_file) as fh:

            for line in fh:
                if re.match(".*Creating operator.*", line):
                    # Retrieves the process name from the string
                    match = re.match(".*Creating operator > (.*) --", line)
                    process = match.group(1)

                    if any([process.startswith(x) for x in self._blacklist]):
                        continue

                    if process not in self.skip_processes:
                        self.processes[match.group(1)] = {
                            "barrier": "W",
                            "submitted": set(),
                            "finished": set(),
                            "failed": set(),
                            "retry": set(),
                            "cpus": None,
                            "memory": None
                        }
                        self.process_tags[process] = {}

                # Retrieves the pipeline name from the string
                if re.match(".*Launching `.*` \[.*\] ", line):
                    tag_match = re.match(".*Launching `.*` \[(.*)\] ", line)
                    self.pipeline_tag = tag_match.group(1) if tag_match else \
                        "?"
                    name_match = re.match(".*Launching `(.*)` \[.*\] ", line)
                    self.pipeline_name = name_match.group(1) if name_match \
                        else "?"

        self.content_lines = len(self.processes)

    def _clear_inspect(self):
        """Clears inspect attributes when re-executing a pipeline"""

        self.trace_info = defaultdict(list)
        self.process_tags = {}
        self.process_stats = {}
        self.samples = []
        self.stored_ids = []
        self.stored_log_ids = []
        self.time_start = None
        self.time_stop = None
        self.execution_command = None
        self.nextflow_version = None
        self.abort_cause = None
        self._c = 0
        # Clean up of tag running status
        for p in self.processes.values():
            p["barrier"] = "W"
            for i in ["submitted", "finished", "failed", "retry"]:
                p[i] = set()

    def _update_pipeline_status(self):
        """Parses the .nextflow.log file for signatures of pipeline status.
        It sets the :attr:`status_info` attribute.
        """

        with open(self.log_file) as fh:

            first_line = next(fh)
            time_str = " ".join(first_line.split()[:2])
            self.time_start = time_str

            if not self.execution_command:
                try:
                    self.execution_command = re.match(
                        ".*nextflow run (.*)", first_line).group(1)
                except AttributeError:
                    self.execution_command = "Unknown"

            for line in fh:

                if "DEBUG nextflow.cli.CmdRun" in line:
                    if not self.nextflow_version:
                        try:
                            vline = next(fh)
                            self.nextflow_version = re.match(
                                ".*Version: (.*)", vline).group(1)
                        except AttributeError:
                            self.nextflow_version = "Unknown"

                if "Session aborted" in line:
                    self.run_status = "aborted"
                    # Get abort cause
                    try:
                        self.abort_cause = re.match(
                            ".*Cause: (.*)", line).group(1)
                    except AttributeError:
                        self.abort_cause = "Unknown"
                    # Get time of pipeline stop
                    time_str = " ".join(line.split()[:2])
                    self.time_stop = time_str
                    self.send = True
                    return
                if "Execution complete -- Goodbye" in line:
                    self.run_status = "complete"
                    # Get time of pipeline stop
                    time_str = " ".join(line.split()[:2])
                    self.time_stop = time_str
                    self.send = True
                    return

        if self.run_status not in ["running", ""]:
            self._clear_inspect()
            # Take a break to allow nextflow to restart before refreshing
            # pipeine processes
            sleep(5)
            self._get_pipeline_processes()

        self.run_status = "running"

    def _update_tag_status(self, process, vals):
        """ Updates the 'submitted', 'finished', 'failed' and 'retry' status
        of each process/tag combination.

        Process/tag combinations provided to this method already appear on
        the trace file, so their submission status is updated based on their
        execution status from nextflow.

        For instance, if a tag is successfully
        complete, it is moved from the 'submitted' to the 'finished' list.
        If not, it is moved to the 'failed' list.

        Parameters
        ----------
        process : str
            Name of the current process. Must be present in attr:`processes`
        vals : list
            List of tags for this process that have been gathered in the
            trace file.
        """

        good_status = ["COMPLETED", "CACHED"]

        # Update status of each process
        for v in list(vals)[::-1]:
            p = self.processes[process]
            tag = v["tag"]

            # If the process/tag is in the submitted list, move it to the
            # complete or failed list
            if tag in p["submitted"]:
                p["submitted"].remove(tag)
                if v["status"] in good_status:
                    p["finished"].add(tag)
                elif v["status"] == "FAILED":
                    if not v["work_dir"]:
                        v["work_dir"] = ""
                    self.process_tags[process][tag]["log"] = \
                        self._retrieve_log(join(v["work_dir"], ".command.log"))
                    p["failed"].add(tag)

            # It the process/tag is in the retry list and it completed
            # successfully, remove it from the retry and fail lists. Otherwise
            # maintain it in the retry/failed lists
            elif tag in p["retry"]:
                if v["status"] in good_status:
                    p["retry"].remove(tag)
                    p["failed"].remove(tag)
                    del self.process_tags[process][tag]["log"]
                elif self.run_status == "aborted":
                    p["retry"].remove(tag)

            elif v["status"] in good_status:
                p["finished"].add(tag)

            # Filter tags without a successfull status.
            if v["status"] not in good_status:
                if v["tag"] in list(p["submitted"]) + list(p["finished"]):
                    vals.remove(v)

        return vals

    def _update_barrier_status(self):
        """Checks whether the channels to each process have been closed.
        """

        with open(self.log_file) as fh:

            for line in fh:

                # Exit barrier update after session abort signal
                if "Session aborted" in line:
                    return

                if "<<< barrier arrive" in line:
                    # Retrieve process name from string
                    process_m = re.match(".*process: (.*)\)", line)
                    if process_m:
                        process = process_m.group(1)
                        # Updates process channel to complete
                        if process in self.processes:
                            self.processes[process]["barrier"] = "C"

    @staticmethod
    def _retrieve_log(path):
        """Method used to retrieve the contents of a log file into a list.

        Parameters
        ----------
        path

        Returns
        -------
        list or None
            Contents of the provided file, each line as a list entry
        """

        if not os.path.exists(path):
            return None

        with open(path) as fh:
            return fh.readlines()

    def _update_trace_info(self, fields, hm):
        """Parses a trace line and updates the :attr:`status_info` attribute.

        Parameters
        ----------
        fields : list
            List of the tab-seperated elements of the trace line
        hm : dict
            Maps the column IDs to their position in the fields argument.
            This dictionary object is retrieve from :func:`_header_mapping`.
        """

        process = fields[hm["process"]]

        if process not in self.processes:
            return

        # Get information from a single line of trace file
        info = dict((column, fields[pos]) for column, pos in hm.items())

        # The headers that will be used to populate the process
        process_tag_headers = ["realtime", "rss", "rchar", "wchar"]
        for h in process_tag_headers:

            # In the rare occasion the tag is parsed first in the trace
            # file than the log file, add the new tag.
            if info["tag"] not in self.process_tags[process]:
                # If the 'start' tag is present in the trace, use that
                # information. If not, it will be parsed in the log file.
                try:
                    timestart = info["start"].split()[1]
                except KeyError:
                    timestart = None
                self.process_tags[process][info["tag"]] = {
                    "workdir": self._expand_path(info["hash"]),
                    "start": timestart
                }

            if h in info and info["tag"] != "-":
                if h != "realtime" and info[h] != "-":
                    self.process_tags[process][info["tag"]][h] = \
                        round(self._size_coverter(info[h]), 2)
                else:
                    self.process_tags[process][info["tag"]][h] = info[h]

        # Set allocated cpu and memory information to process
        if "cpus" in info and not self.processes[process]["cpus"]:
            self.processes[process]["cpus"] = info["cpus"]
        if "memory" in info and not self.processes[process]["memory"]:
            try:
                self.processes[process]["memory"] = self._size_coverter(
                    info["memory"])
            except ValueError:
                self.processes[process]["memory"] = None

        if info["hash"] in self.stored_ids:
            return

        # If the task hash code is provided, expand it to the work directory
        # and add a new entry
        if "hash" in info:
            hs = info["hash"]
            info["work_dir"] = self._expand_path(hs)

        if "tag" in info:
            tag = info["tag"]
            if tag != "-" and tag not in self.samples and \
                    tag.split()[0] not in self.samples:
                self.samples.append(tag)

        self.trace_info[process].append(info)
        self.stored_ids.append(info["hash"])

    def _update_process_resources(self, process, vals):
        """Updates the resources info in :attr:`processes` dictionary.
        """

        resources = ["cpus"]

        for r in resources:
            if not self.processes[process][r]:
                try:
                    self.processes[process][r] = vals[0]["cpus"]
                # When the trace column is not present
                except KeyError:
                    pass

    def _cpu_load_parser(self, cpus, cpu_per, t):
        """Parses the cpu load from the number of cpus and its usage
        percentage and returnsde cpu/hour measure

        Parameters
        ----------
        cpus : str
            Number of cpus allocated.
        cpu_per : str
            Percentage of cpu load measured (e.g.: 200,5%).
        t : str
            The time string can be something like '20s', '1m30s' or '300ms'.
        """

        try:
            _cpus = float(cpus)
            _cpu_per = float(cpu_per.replace(",", ".").replace("%", ""))
            hours = self._hms(t) / 60 / 24

            return ((_cpu_per / (100 * _cpus)) * _cpus) * hours

        except ValueError:
            return 0

    def _assess_resource_warnings(self, process, vals):
        """Assess whether the cpu load or memory usage is above the allocation

        Parameters
        ----------
        process : str
            Process name
        vals : vals
            List of trace information for each tag of that process

        Returns
        -------
        cpu_warnings : dict
            Keys are tags and values are the excessive cpu load
        mem_warnings : dict
            Keys are tags and values are the excessive rss
        """

        cpu_warnings = {}
        mem_warnings = {}

        for i in vals:
            try:
                expected_load = float(i["cpus"]) * 100
                cpu_load = float(i["%cpu"].replace(",", ".").replace("%", ""))

                if expected_load * 0.9 > cpu_load > expected_load * 1.10:
                    cpu_warnings[i["tag"]] = {
                        "expected":  expected_load,
                        "value": cpu_load
                    }
            except (ValueError, KeyError):
                pass

            try:
                rss = self._size_coverter(i["rss"])
                mem_allocated = self._size_coverter(i["memory"])

                if rss > mem_allocated * 1.10:
                    mem_warnings[i["tag"]] = {
                        "expected": mem_allocated,
                        "value": rss
                    }
            except (ValueError, KeyError):
                pass

        return cpu_warnings, mem_warnings

    def _update_process_stats(self):
        """Updates the process stats with the information from the processes

        This method is called at the end of each static parsing of the nextflow
        trace file. It re-populates the :attr:`process_stats` dictionary
        with the new stat metrics.
        """

        good_status = ["COMPLETED", "CACHED"]

        for process, vals in self.trace_info.items():

            # Update submission status of tags for each process
            vals = self._update_tag_status(process, vals)

            # Update process resources
            self._update_process_resources(process, vals)

            self.process_stats[process] = {}

            inst = self.process_stats[process]

            # Get number of completed samples
            inst["completed"] = "{}".format(
                len([x for x in vals if x["status"] in good_status]))

            # Get average time
            try:
                time_array = [self._hms(x["realtime"]) for x in vals]
                mean_time = round(sum(time_array) / len(time_array), 1)
                mean_time_str = strftime('%H:%M:%S', gmtime(mean_time))
                inst["realtime"] = mean_time_str
            # When the realtime column is not present
            except KeyError:
                inst["realtime"] = "-"

            # Get cumulative cpu/hours
            try:
                cpu_hours = [self._cpu_load_parser(
                    x["cpus"], x["%cpu"], x["realtime"]) for x in vals]
                inst["cpuhour"] = round(sum(cpu_hours), 2)
            # When the realtime, cpus or %cpus column are not present
            except KeyError:
                inst["cpuhour"] = "-"

            # Assess resource warnings
            inst["cpu_warnings"], inst["mem_warnings"] = \
                self._assess_resource_warnings(process, vals)

            # Get maximum memory
            try:
                rss_values = [self._size_coverter(x["rss"]) for x in vals
                              if x["rss"] != "-"]
                if rss_values:
                    max_rss = round(max(rss_values))
                    rss_str = self._size_compress(max_rss)
                else:
                    rss_str = "-"
                inst["maxmem"] = rss_str
            except KeyError:
                inst["maxmem"] = "-"

            # Get read size
            try:
                rchar_values = [self._size_coverter(x["rchar"]) for x in vals
                                if x["rchar"] != "-"]
                if rchar_values:
                    avg_rchar = round(sum(rchar_values) / len(rchar_values))
                    rchar_str = self._size_compress(avg_rchar)
                else:
                    rchar_str = "-"
            except KeyError:
                rchar_str = "-"
            inst["avgread"] = rchar_str

            # Get write size
            try:
                wchar_values = [self._size_coverter(x["wchar"]) for x in vals
                                if x["wchar"] != "-"]
                if wchar_values:
                    avg_wchar = round(sum(wchar_values) / len(wchar_values))
                    wchar_str = self._size_compress(avg_wchar)
                else:
                    wchar_str = "-"
            except KeyError:
                wchar_str = "-"
            inst["avgwrite"] = wchar_str

    #################
    # PARSING METHODS
    #################

    def trace_parser(self):
        """Method that parses the trace file once and updates the
        :attr:`status_info` attribute with the new entries.
        """

        # Check the timestamp of the tracefile. Only proceed with the parsing
        # if it changed from the previous time.
        size_stamp = os.path.getsize(self.trace_file)
        self.trace_retry = 0
        if size_stamp and size_stamp == self.trace_sizestamp:
            return
        else:
            logger.debug("Updating trace size stamp to: {}".format(size_stamp))
            self.trace_sizestamp = size_stamp

        with open(self.trace_file) as fh:

            # Skip potential empty lines at the start of file
            header = next(fh).strip()
            while not header:
                header = next(fh).strip()

            # Get header mappings before parsing the file
            hm = self._header_mapping(header)

            for line in fh:

                # Skip empty lines
                if line.strip() == "":
                    continue

                fields = line.strip().split("\t")

                # Skip if task ID was already processes
                if fields[hm["task_id"]] in self.stored_ids:
                    continue

                # Parse trace entry and update status_info attribute
                self._update_trace_info(fields, hm)
                self.send = True

        self._update_process_stats()
        self._update_barrier_status()

    def log_parser(self):
        """Method that parses the nextflow log file once and updates the
        submitted number of samples for each process
        """

        # Check the timestamp of the log file. Only proceed with the parsing
        # if it changed from the previous time.
        size_stamp = os.path.getsize(self.log_file)
        self.log_retry = 0
        if size_stamp and size_stamp == self.log_sizestamp:
            return
        else:
            logger.debug("Updating log size stamp to: {}".format(size_stamp))
            self.log_sizestamp = size_stamp

        # Regular expression to catch four groups:
        # 1. Start timestamp
        # 2. Work directory hash
        # 3. Process name
        # 4. Tag name
        r = ".* (.*) \[.*\].*\[(.*)\].*process > (.*) \((.*)\).*"

        with open(self.log_file) as fh:

            for line in fh:
                if "Submitted process >" in line or \
                        "Re-submitted process >" in line or \
                        "Cached process >" in line:
                    m = re.match(r, line)
                    if not m:
                        continue

                    time_start = m.group(1)
                    workdir = m.group(2)
                    process = m.group(3)
                    tag = m.group(4)

                    # Skip if this line has already been parsed
                    if time_start + tag not in self.stored_log_ids:
                        self.stored_log_ids.append(time_start + tag)
                    else:
                        continue

                    # For first time processes
                    if process not in self.processes:
                        continue
                    p = self.processes[process]

                    # Skip is process/tag combination has finished or is retrying
                    if tag in list(p["finished"]) + list(p["retry"]):
                        continue

                    # Update failed process/tags when they have been re-submitted
                    if tag in list(p["failed"]) and \
                            "Re-submitted process >" in line:
                        p["retry"].add(tag)
                        self.send = True
                        continue

                    # Set process barrier to running. Check for barrier status
                    # are performed at the end of the trace parsing in the
                    # _update_barrier_status method.
                    p["barrier"] = "R"
                    if tag not in p["submitted"]:
                        p["submitted"].add(tag)
                        # Update the process_tags attribute with the new tag.
                        # Update only when the tag does not exist. This may rarely
                        # occur when the tag is parsed first in the trace file
                        if tag not in self.process_tags[process]:
                            self.process_tags[process][tag] = {
                                "workdir": self._expand_path(workdir),
                                "start": time_start
                            }
                            self.send = True
                        # When the tag is filled in the trace file parsing,
                        # the timestamp may not be present in the trace. In
                        # those cases, fill that information here.
                        elif not self.process_tags[process][tag]["start"]:
                            self.process_tags[process][tag]["start"] = time_start
                            self.send = True

        self._update_pipeline_status()

    def update_inspection(self):
        """Wrapper method that calls the appropriate main updating methods of
        the inspection.

        It is meant to be used inside a loop (like while), so that it can
        continuously update the class attributes from the trace and log files.
        It already implements checks to parse these files only when they
        change, and they ignore entries that have been previously processes.
        """

        try:
            self.log_parser()
        except (FileNotFoundError, StopIteration) as e:
            logger.debug("ERROR: " + str(sys.exc_info()[0]))
            self.log_retry += 1
            if self.log_retry == self.MAX_RETRIES:
                raise e
        try:
            self.trace_parser()
        except (FileNotFoundError, StopIteration) as e:
            logger.debug("ERROR: " + str(sys.exc_info()[0]))
            self.trace_retry += 1
            if self.trace_retry == self.MAX_RETRIES:
                raise e

    #################
    # CURSES METHODS
    #################

    def display_overview(self):
        """Displays the default pipeline inspection overview
        """

        stay_alive = True

        self.screen = curses.initscr()

        self.screen.keypad(True)
        self.screen.nodelay(-1)
        curses.cbreak()
        curses.noecho()
        curses.start_color()

        self.screen_lines = self.screen.getmaxyx()[0]
        # self.screen_width = self.screen.getmaxyx()[1]

        try:
            while stay_alive:

                # Provide functionality to certain keybindings
                self._curses_keybindings()
                # Updates main inspector attributes
                self.update_inspection()
                # Display curses interface
                self.flush_overview()

                sleep(self.refresh_rate)
        except FileNotFoundError:
            sys.stderr.write(colored_print(
                "ERROR: nextflow log and/or trace files are no longer "
                "reachable!", "red_bold"))
        except Exception as e:
            sys.stderr.write(str(e))
        finally:
            curses.nocbreak()
            self.screen.keypad(0)
            curses.echo()
            curses.endwin()

    def _curses_keybindings(self):

        c = self.screen.getch()
        # Provide scroll up/down with keys or mouse wheel
        if c == curses.KEY_UP:
            self._updown("up")
        elif c == curses.KEY_DOWN:
            self._updown("down")
        elif c == curses.KEY_LEFT:
            self._rightleft("left")
        elif c == curses.KEY_RIGHT:
            self._rightleft("right")
        # Trigger screen size update on resize
        elif c == curses.KEY_RESIZE:
            self.screen_lines = self.screen.getmaxyx()[0]
        # Exit interface when pressing q
        elif c == ord('q'):
            raise Exception

    def _updown(self, direction):
        """Provides curses scroll functionality.
        """

        if direction == "up" and self.top_line != 0:
            self.top_line -= 1
        elif direction == "down" and \
                self.screen.getmaxyx()[0] + self.top_line\
                <= self.content_lines + 3:
            self.top_line += 1

    def _rightleft(self, direction):
        """Provides curses horizontal padding"""

        if direction == "left" and self.padding != 0:
            self.padding -= 1

        if direction == "right" and \
                self.screen.getmaxyx()[1] + self.padding < self.max_width:
            self.padding += 1

    def flush_overview(self):
        """Displays the default overview of the pipeline execution from the
        :attr:`status_info`, :attr:`processes` and :attr:`run_status`
        attributes into stdout.
        """

        colors = {
            "W": 1,
            "R": 2,
            "C": 3
        }

        pc = {
            "running": 3,
            "complete": 3,
            "aborted": 4,
            "error": 4
        }

        curses.init_pair(1, curses.COLOR_WHITE, curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_BLUE, curses.COLOR_BLACK)
        curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)
        curses.init_pair(4, curses.COLOR_MAGENTA, curses.COLOR_BLACK)

        # self.screen.erase()

        height, width = self.screen.getmaxyx()
        win = curses.newpad(height, 2000)

        # Add static header
        header = "Pipeline [{}] inspection at {}. Status: ".format(
            self.pipeline_tag, strftime("%Y-%m-%d %H:%M:%S", gmtime()))

        win.addstr(0, 0, header)
        win.addstr(0, len(header), self.run_status,
                   curses.color_pair(pc[self.run_status]))
        submission_str = "{0:23.23}  {1:23.23}  {2:23.23}  {3:23.23}".format(
            "Running: {}".format(
                sum([len(x["submitted"]) for x in self.processes.values()])
            ),
            "Failed: {}".format(
                sum([len(x["failed"]) for x in self.processes.values()])
            ),
            "Retrying: {}".format(
                sum([len(x["retry"]) for x in self.processes.values()])
            ),
            "Completed: {}".format(
                sum([len(x["finished"]) for x in self.processes.values()])
            )
        )

        win.addstr(
            1, 0, submission_str, curses.color_pair(1)
        )

        headers = ["", "Process", "Running", "Complete", "Error",
                   "Avg Time", "Max Mem", "Avg Read", "Avg Write"]
        header_str = "{0: ^1} " \
                     "{1: ^25}  " \
                     "{2: ^7} " \
                     "{3: ^7} " \
                     "{4: ^7} " \
                     "{5: ^10} " \
                     "{6: ^10} " \
                     "{7: ^10} " \
                     "{8: ^10} ".format(*headers)
        self.max_width = len(header_str)
        win.addstr(3, 0, header_str, curses.A_UNDERLINE | curses.A_REVERSE)

        # Get display size
        top = self.top_line
        bottom = self.screen_lines - 4 + self.top_line

        # Fetch process information
        for p, process in enumerate(
                list(self.processes.keys())[top:bottom]):

            if process not in self.process_stats:
                vals = ["-"] * 8
                txt_fmt = curses.A_NORMAL
            else:
                ref = self.process_stats[process]
                vals = [ref["completed"],
                        len(self.processes[process]["failed"]),
                        ref["realtime"],
                        ref["maxmem"], ref["avgread"],
                        ref["avgwrite"]]
                txt_fmt = curses.A_BOLD

            proc = self.processes[process]
            if proc["retry"]:
                completed = "{}({})".format(len(proc["submitted"]),
                                            len(proc["retry"]))
            else:
                completed = "{}".format(len(proc["submitted"]))

            win.addstr(
                4 + p, 0, "{0: ^1} "
                          "{1:25.25}  "
                          "{2: ^7} "
                          "{3: ^7} "
                          "{4: ^7} "
                          "{5: ^10} "
                          "{6: ^10} "
                          "{7: ^10} "
                          "{8: ^10} ".format(
                                proc["barrier"],
                                process,
                                completed,
                                *vals),
                curses.color_pair(colors[proc["barrier"]]) | txt_fmt)

        win.clrtoeol()
        win.refresh(0, self.padding, 0, 0, height-1, width-1)

    ###################
    # BROADCAST METHODS
    ###################

    def _convert_process_dict(self):

        d = {}

        for k, v in self.processes.items():
            d[k] = {
                "barrier": v["barrier"],
                "cpus": v["cpus"],
                "memory": v["memory"]
            }
            for i in ["submitted", "finished", "failed", "retry"]:
                d[k][i] = list(v[i])

        return d

    def _prepare_table_data(self):

        # Set data mappings
        mappings = {
            "Barrier": "barrier",
            "Process": "process",
            "Running": "running",
            "Complete": "complete",
            "Error": "error",
            "Avg Time": "avgTime",
            "CPU/hour": "cpuhour",
            "Max Mem": "maxMem",
            "Avg Read": "avgRead",
            "Avg Write": "avgWrite"
        }

        # Set table data
        data = []
        table_headers = ["avgTime", "cpuhour", "maxMem", "avgRead", "avgWrite"]
        for process in list(self.processes):

            proc = self.processes[process]
            # Add general data that is always available for all processes
            current_data = {
                "process": process,
                "barrier": proc["barrier"],
                "complete": list(proc["finished"]),
                "error": list(proc["failed"]),
                "running": list(proc["submitted"])
            }

            # Add stats data that is only available for processes that have
            # finished once.
            if process not in self.process_stats:
                current_data = {
                    **current_data,
                    **dict((x, "-") for x in table_headers),
                    **{"cpuWarn": {}, "memWarn": {}}
                }

            else:
                ref = self.process_stats[process]
                current_data = {
                    **current_data,
                    **{"avgTime": ref["realtime"],
                       "cpuhour": ref["cpuhour"],
                       "maxMem": ref["maxmem"],
                       "avgRead": ref["avgread"],
                       "avgWrite": ref["avgwrite"],
                       "cpuWarn": ref["cpu_warnings"],
                       "memWarn": ref["mem_warnings"]}
                }

            data.append(current_data)

        return mappings, data

    def _prepare_overview_data(self):

        return [
            {
                "header": "Pipeline name",
                "value": self.pipeline_name
            },
            {
                "header": "Pipeline tag",
                "value": self.pipeline_tag
            },
            {
                "header": "Number of processes",
                "value": len(self.processes)
            }]

    def _prepare_general_details(self):
        return [
            {
                "header": "Pipeline directory",
                "value": self.workdir
            },
            {
                "header": "Work directory",
                "value": join(self.workdir, "work")
            },
            {
                "header": "Nextflow command",
                "value": self.execution_command
            },
            {
                "header": "Nextflow version",
                "value": self.nextflow_version
            }
        ]

    def _get_log_lines(self, n=300):
        """Returns a list with the last ``n`` lines of the nextflow log file

        Parameters
        ----------
        n : int
            Number of last lines from the log file

        Returns
        -------
        list
            List of strings with the nextflow log
        """

        with open(self.log_file) as fh:
            last_lines = fh.readlines()[-n:]

        return last_lines

    def _prepare_run_status_data(self):

        if self.run_status == "aborted":
            log_lines = self._get_log_lines()
        else:
            log_lines = None

        return {
            "value": self.run_status,
            "abortCause": self.abort_cause,
            "logLines": log_lines
        }

    def _send_status_info(self, run_id):

        mappings, data = self._prepare_table_data()
        overview_data = self._prepare_overview_data()
        general_details = self._prepare_general_details()
        status_data = self._prepare_run_status_data()

        status_json = {
            "generalOverview": overview_data,
            "generalDetails": general_details,
            "tableData": data,
            "tableMappings": mappings,
            "processInfo": self._convert_process_dict(),
            "processTags": self.process_tags,
            "runStatus": status_data,
            "timeStart": str(self.time_start),
            "timeStop": str(self.time_stop) if self.time_stop else "-",
            "processes": list(self.processes)
        }

        self._c += 1
        logger.debug("Payload [{}] sent with size: {}".format(
            self._c,
            asizeof.asizeof(json.dumps(status_json))
        ))

        try:
            requests.put(self.broadcast_address,
                         json={"run_id": run_id, "status_json": status_json})
        except requests.exceptions.ConnectionError:
            logger.error(colored_print(
                "ERROR: Could not establish connection with server. The server"
                " may be down or there is a problem with your internet "
                "connection.", "red_bold"))
            sys.exit(1)

    def _prepare_static_info(self):
        """Prepares the first batch of information, containing static
        information such as the pipeline file, and configuration files

        Returns
        -------
        dict
            Dict with the static information for the first POST request
        """

        pipeline_files = {}

        with open(join(self.workdir, self.pipeline_name)) as fh:
            pipeline_files["pipelineFile"] = fh.readlines()

        nf_config = join(self.workdir, "nextflow.config")
        if os.path.exists(nf_config):
            with open(nf_config) as fh:
                pipeline_files["configFile"] = fh.readlines()

        # Check for specific flowcraft configurations files
        configs = {
            "params.config": "paramsFile",
            "resources.config": "resourcesFile",
            "containers.config": "containersFile",
            "user.config": "userFile",
        }
        for config, key in configs.items():
            cfile = join(self.workdir, config)
            if os.path.exists(cfile):
                with open(cfile) as fh:
                    pipeline_files[key] = fh.readlines()

        return pipeline_files

    def _dag_file_to_dict(self):
        """Function that opens the dotfile named .treeDag.json in the current
        working directory

        Returns
        -------
        Returns a dictionary with the dag object to be used in the post
        instance available through the method _establish_connection

        """
        try:
            dag_file = open(os.path.join(self.workdir, ".treeDag.json"))
            dag_json = json.load(dag_file)
        except (FileNotFoundError, json.decoder.JSONDecodeError):
            logger.warning(colored_print(
                "WARNING: dotfile named .treeDag.json not found or corrupted",
                "red_bold"))
            dag_json = {}

        return dag_json

    def _establish_connection(self, run_id, dict_dag):

        try:

            static_info = self._prepare_static_info()

            logger.debug("Sending initial data with run id: {}".format(run_id))

            payload = {"run_id": run_id, "dag_json": dict_dag,
                       "pipeline_files": static_info}
            logger.debug("Connection payload size: {}".format(
                asizeof.asizeof(payload)))

            r = requests.post(self.broadcast_address,
                              json=payload)

            logger.debug("Response received: {}".format(r.status_code))
            if r.status_code != 201:
                logger.error(colored_print(
                    "ERROR: There was a problem sending data to the server"
                    "with reason: {}".format(r.reason)))
                sys.exit(1)
        except requests.exceptions.ConnectionError:
            logger.error(colored_print(
                "ERROR: Could not establish connection with server. The server"
                " may be down or there is a problem with your internet "
                "connection.", "red_bold"))
            sys.exit(1)

    def _close_connection(self, run_id):

        try:
            r = requests.delete(self.broadcast_address,
                                json={"run_id": run_id})
            if r.status_code != 202:
                logger.error(colored_print(
                    "ERROR: There was a problem sending data to the server"
                    "with reason: {}".format(r.reason)))
        except requests.exceptions.ConnectionError:
            logger.error(colored_print(
                "ERROR: Could not establish connection with server. The server"
                " may be down or there is a problem with your internet "
                "connection.", "red_bold"))
            sys.exit(1)

    def _get_run_hash(self):
        """Gets the hash of the nextflow file"""

        # Get name of the pipeline from the log file
        with open(self.log_file) as fh:
            header = fh.readline()
        pipeline_path = re.match(".*nextflow run ([^\s]+).*", header).group(1)

        # Get hash from the entire pipeline file
        pipeline_hash = hashlib.md5()
        with open(pipeline_path, "rb") as fh:
            for chunk in iter(lambda: fh.read(4096), b""):
                pipeline_hash.update(chunk)
        # Get hash from the current working dir and hostname
        workdir = self.workdir.encode("utf8")
        hostname = socket.gethostname().encode("utf8")
        dir_hash = hashlib.md5(workdir + hostname)

        return pipeline_hash.hexdigest() + dir_hash.hexdigest()

    def _print_msg(self, run_id):

        inspect_address = "{}inspect/{}".format(self.app_address, run_id)
        logger.info(colored_print(
            "Starting broadcast. You can see the pipeline progress on the "
            "link below:", "green_bold"))
        logger.info("{}".format(inspect_address))

    def broadcast_status(self):

        logger.info(colored_print("Preparing broadcast data...", "green_bold"))

        run_hash = self._get_run_hash()
        dict_dag = self._dag_file_to_dict()
        _broadcast_sent = False
        logger.debug("Establishing connection...")
        self._establish_connection(run_hash, dict_dag)

        stay_alive = True
        try:
            logger.debug("Starting inspection loop")
            while stay_alive:

                if not _broadcast_sent:
                    self._print_msg(run_hash)
                    _broadcast_sent = True

                self.update_inspection()
                if self.send:
                    logger.debug("Updating inspection")
                    self._send_status_info(run_hash)
                    self.send = False

                sleep(self.refresh_rate)

        except FileNotFoundError:
            logger.error(colored_print(
                "ERROR: nextflow log and/or trace files are no longer "
                "reachable!", "red_bold"))
        except Exception:
            logger.exception("ERROR: " + str(sys.exc_info()[0]))
        finally:
            logger.info("Closing connection")
            self._close_connection(run_hash)
