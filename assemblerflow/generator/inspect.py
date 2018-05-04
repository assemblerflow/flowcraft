import re
import os
import sys
import curses
import signal
import locale

from os.path import join, abspath
from time import gmtime, strftime, sleep
from collections import defaultdict, OrderedDict

try:
    import generator.error_handling as eh
except ImportError:
    import assemblerflow.generator.error_handling as eh

locale.setlocale(locale.LC_ALL, '')
code = locale.getpreferredencoding()


def signal_handler(screen):
    """This function is bound to the SIGINT signal (like ctrl+c) to graciously
    exit the program and reset the curses options.
    """

    screen.clear()
    screen.refresh()

    curses.nocbreak()
    screen.keypad(0)
    curses.echo()
    curses.endwin()
    print("Exiting assemblerflow inspection... Bye")
    sys.exit(0)


class NextflowInspector:

    def __init__(self, trace_file, refresh_rate, pretty=False):

        self.trace_file = trace_file
        """
        str: Path to nextflow trace file.
        """

        self.trace_sizestamp = None
        """
        str: Stores the timestamp of the last modification of the trace file.
        This is used to parse the file only when it has changed.
        """

        self.refresh_rate = refresh_rate
        """
        float: Frequency (in seconds) that the curses screen will be refreshed.
        """

        self.stored_ids = []
        """
        list: Stores the task_ids that have already been parsed.
        """

        self.process_info = defaultdict(list)
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
        :func:`_update_process_stats`.
        """

        self.samples = []
        """
        list: List of samples inferred from the pipeline.
        """

        self.skip_processes = ["status", "compile_status", "report",
                               "compile_reports", "fullConsensus"]
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
        str: Stores the timestamp of the last modification of the nextflow 
        log file. This is used to parse the file only when it has changed.
        """

        self.pipeline_name = ""
        """
        str: Name of the pipeline, parsed from .nextflow.log
        """

        self.run_status = ""
        """
        str: Status of the pipeline. Can be either 'running', 'aborted',
        'error', 'complete'.
        """

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
        self.screen_lines = None
        self.content_lines = 0

        self._check_required_files()
        self._parser_pipeline_processes()
        self._get_pipeline_status()

        # Bind SIGINT to singal_handler function
        signal.signal(signal.SIGINT, lambda *x: signal_handler(self.screen))

    #################
    # UTILITY METHODS
    #################

    def _check_required_files(self):

        if not os.path.exists(self.trace_file):
            raise eh.InspectionError("The provided trace file could not be "
                                     "opened: {}".format(self.trace_file))

        if not os.path.exists(self.log_file):
            raise eh.InspectionError("The .nextflow.log files could not be "
                                     "opened. Are you sure you are in a "
                                     "nextflow project directory?")

    def _get_pipeline_status(self):
        """Parses the .nextflow.log file for signatures of pipeline status.
        It sets the :attr:`status_info` attribute.
        """

        with open(self.log_file) as fh:

            for line in fh:
                if "Session aborted" in line:
                    self.run_status = "aborted"
                    return
                if "Execution complete -- Goodbye" in line:
                    self.run_status = "complete"
                    return

        self.run_status = "running"

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

        return dict((x, pos) for pos, x in enumerate(header.split("\t")))

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

        first_hash, second_hash = hash_str.split("/")
        first_hash_path = join(abspath("work"), first_hash)

        for l in os.listdir(first_hash_path):
            if l.startswith(second_hash):
                return join(first_hash_path, l)

    @staticmethod
    def hms(s):

        if s == "-":
            return 0

        if s.endswith("ms"):
            return int(s.rstrip("ms")) / 1000

        fields = list(map(float, re.split("[hms]", s)[:-1]))
        if len(fields) == 3:
            return fields[0] * 3600 + fields[1] * 60 + fields[2]
        elif len(fields) == 2:
            return fields[0] * 60 + fields[1]
        else:
            return fields[0]

    @staticmethod
    def size_coverter(s):

        if s.endswith("KB"):
            return float(s.rstrip("KB")) / 1024

        elif s.endswith(" B"):
            return float(s.rstrip("B")) / 1024 / 1024

        elif s.endswith("MB"):
            return float(s.rstrip("MB"))

        elif s.endswith("GB"):
            return float(s.rstrip("GB")) * 1024

        else:
            return float(s)

    @staticmethod
    def size_compress(s):

        if s / 1024 > 1:
            return "{}GB".format(round(s / 1024), 1)
        else:
            return "{}MB".format(s)

    #################
    # PARSE METHODS
    #################

    def _parser_pipeline_processes(self):
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
                    match = re.match(".*Creating operator > (.*) --", line)
                    process = match.group(1)

                    if any([process.startswith(x) for x in self._blacklist]):
                        continue

                    if process not in self.skip_processes:
                        self.processes[match.group(1)] = {
                            "barrier": "W",
                            "submitted": [],
                            "finished": []
                        }

                if re.match(".*Launching `.*` \[.*\] ", line):
                    match = re.match(".*Launching `.*` \[(.*)\] ", line)
                    self.pipeline_name = match.group(1)

        self.content_lines = len(self.processes)

    def _update_barrier_status(self):
        """Updates the run_status key of the :attr:`process_stats` from the
        config.
        """

        with open(self.log_file) as fh:

            for line in fh:

                # Exit barrier update after session abort signal
                if "Session aborted" in line:
                    return

                if "<<< barrier arrive" in line:
                    process_m = re.match(".*process: (.*)\)", line)
                    if process_m:
                        process = process_m.group(1)

                        if process in self.processes:
                            self.processes[process]["barrier"] = "C"

    def _update_status(self, fields, hm):
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

        self.process_info[process].append(info)
        self.stored_ids.append(info["hash"])

    def static_parser(self):
        """Method that parses the trace file once and updates the
        :attr:`status_info` attribute with the new entries.
        """

        # Check the timestamp of the tracefile. Only proceed with the parsing
        # if it changed from the previous time.
        size_stamp = os.path.getsize(self.trace_file)
        if size_stamp and size_stamp == self.trace_sizestamp:
            return
        else:
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
                self._update_status(fields, hm)

        self._update_process_stats()
        self._get_pipeline_status()
        self._update_barrier_status()

    def static_log_parser(self):
        """Method that parses the nextflow log file once and updates the
        submitted number of samples for each process
        """

        # Check the timestamp of the log file. Only proceed with the parsing
        # if it changed from the previous time.
        size_stamp = os.path.getsize(self.log_file)
        if size_stamp and size_stamp == self.log_sizestamp:
            return
        else:
            self.log_sizestamp = size_stamp

        with open(self.log_file) as fh:

            for line in fh:
                if "Submitted process >" in line or \
                        "Re-submitted process >" in line:
                    m = re.match(".*process > (.*) \((.*)\).*", line)
                    if not m:
                        continue

                    process = m.group(1)
                    sample = m.group(2)

                    if process not in self.processes:
                        continue
                    if sample in self.processes[process]["finished"]:
                        continue
                    if sample in self.processes[process]["submitted"]:
                        continue

                    self.processes[process]["barrier"] = "R"
                    self.processes[process]["submitted"].append(sample)

    def _update_process_stats(self):
        """Updates the process stats with the information from the processes

        This method is called at the end of each static parsing of the nextflow
        trace file. It re-populates the :attr:`process_stats` dictionary
        with the new stat metrics.
        """

        good_status = ["COMPLETED", "CACHED"]

        for process, vals in self.process_info.items():

            self.process_stats[process] = {}

            inst = self.process_stats[process]

            # Remove from the submitted samples
            for v in list(vals):
                p = self.processes[process]
                if v["tag"] in p["submitted"]:
                    p["submitted"].remove(v["tag"])
                    p["finished"].append(v["tag"])
                if v["status"] not in good_status:
                    if v["tag"] in p["submitted"] + p["finished"]:
                        vals.remove(v)

            # Get number of completed samples
            inst["completed"] = "{}".format(
                len([x for x in vals if x["status"] in good_status]))

            # Get number of bad samples
            inst["bad_samples"] = "{}".format(
                len([x for x in vals if x["status"] not in good_status]))

            # Get average time
            time_array = [self.hms(x["realtime"]) for x in vals]
            mean_time = round(sum(time_array) / len(time_array), 1)
            mean_time_str = strftime('%H:%M:%S', gmtime(mean_time))
            inst["realtime"] = mean_time_str

            # Get cumulated time
            # cum_time_str = strftime('%H:%M:%S', gmtime(
            #     round(sum(time_array), 1)))
            # inst["cumtime"] = cum_time_str

            # Get maximum memory
            rss_values = [self.size_coverter(x["rss"]) for x in vals
                          if x["rss"] != "-"]
            if rss_values:
                max_rss = round(max(rss_values))
                rss_str = self.size_compress(max_rss)
            else:
                rss_str = "-"
            inst["maxmem"] = rss_str

            # Get read size
            try:
                rchar_values = [self.size_coverter(x["rchar"]) for x in vals
                                if x["rchar"] != "-"]
                if rchar_values:
                    avg_rchar = round(sum(rchar_values) / len(rchar_values))
                    rchar_str = self.size_compress(avg_rchar)
                else:
                    rchar_str = "-"
            except KeyError:
                rchar_str = "-"
            inst["avgread"] = rchar_str

            # Get write size
            try:
                wchar_values = [self.size_coverter(x["wchar"]) for x in vals
                                if x["wchar"] != "-"]
                if wchar_values:
                    avg_wchar = round(sum(wchar_values) / len(wchar_values))
                    wchar_str = self.size_compress(avg_wchar)
                else:
                    wchar_str = "-"
            except KeyError:
                wchar_str = "-"
            inst["avgwrite"] = wchar_str

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

        try:
            while stay_alive:

                c = self.screen.getch()
                if c == curses.KEY_UP:
                    self.updown("up")
                elif c == curses.KEY_DOWN:
                    self.updown("down")
                elif c == curses.KEY_RESIZE:
                    self.screen_lines = self.screen.getmaxyx()[0]

                self.static_log_parser()
                self.static_parser()
                self.flush_overview()

                sleep(0.1)
        except Exception as e:
            sys.stderr.write(e.msg)
        finally:
            curses.nocbreak()
            self.screen.keypad(0)
            curses.echo()
            curses.endwin()

    def updown(self, direction):
        """Provides curses scroll functionality.
        """

        if direction == "up" and self.top_line != 0:
            self.top_line -= 1
        elif direction == "down" and \
                self.screen.getmaxyx()[0] + self.top_line <= self.content_lines + 5:
            self.top_line += 1

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

        self.screen.erase()

        # Add static header
        header = "Pipeline [{}] inspection at {}. Status: ".format(
            self.pipeline_name, strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        self.screen.addstr(1, 0, header)
        self.screen.addstr(1, len(header), self.run_status,
                           curses.color_pair(pc[self.run_status]))
        self.screen.addstr(
            2, 0, "Running processes: {}".format(
                sum([len(x["submitted"]) for x in self.processes.values()])),
            curses.color_pair(2)
        )
        self.screen.addstr(
            3, 0, "Complete processes: {}".format(
                sum([len(x["finished"]) for x in self.processes.values()])),
            curses.color_pair(3)
        )
        headers = ["", "Process", "Running", "Complete", "Error",
                   "Avg Time", "Max Mem", "Avg Read", "Avg Write"]
        self.screen.addstr(5, 0, "{0: ^1} "
                                 "{1: ^25}  "
                                 "{2: ^7} "
                                 "{3: ^7} "
                                 "{4: ^7} "
                                 "{5: ^10} "
                                 "{6: ^10} "
                                 "{7: ^10} "
                                 "{8: ^10} ".format(*headers),
                           curses.A_UNDERLINE | curses.A_REVERSE)

        # Get display size
        top = self.top_line
        bottom = self.screen_lines - 6 + self.top_line

        # Fetch process information
        for p, process in enumerate(
                list(self.processes.keys())[top:bottom]):

            if process not in self.process_stats:
                vals = ["-"] * 8
                txt_fmt = curses.A_DIM
            else:
                ref = self.process_stats[process]
                vals = [ref["completed"], ref["bad_samples"], ref["realtime"],
                        ref["maxmem"], ref["avgread"],
                        ref["avgwrite"]]
                txt_fmt = curses.A_BOLD

            self.screen.addstr(
                6 + p, 0, "{0: ^1} "
                          "{1:25.25}  "
                          "{2: ^7} "
                          "{3: ^7} "
                          "{4: ^7} "
                          "{5: ^10} "
                          "{6: ^10} "
                          "{7: ^10} "
                          "{8: ^10} ".format(
                                self.processes[process]["barrier"],
                                process,
                                len(self.processes[process]["submitted"]),
                                *vals
            ), curses.color_pair(colors[self.processes[process]["barrier"]]) |
               txt_fmt)

        self.screen.refresh()
