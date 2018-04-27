import re
import os
import sys
import curses
import signal

from os.path import join, abspath
from time import gmtime, strftime, sleep
from collections import defaultdict

# Init curses screen
screen = curses.initscr()

def signal_handler():
    """This function is bound to the SIGINT signal (like ctrl+c) to graciously
    exit the program.
    """

    screen.clear()
    screen.refresh()

    curses.nocbreak()
    screen.keypad(0)
    curses.echo()
    curses.endwin()
    print("Exiting assemblerflow inspection... Bye")
    sys.exit(0)


# Bind SIGINT to singal_handler function
signal.signal(signal.SIGINT, lambda *x: signal_handler())


class NextflowInspector:

    def __init__(self, trace_file, refresh_rate):

        self.trace_file = trace_file
        """
        str: Path to nextflow trace file.
        """

        self.trace_stamp = None
        """
        str: Stores the timestamp of the last modification of the trace file.
        This is used to parse the file only when it has changed.
        """

        self.refresh_rate = refresh_rate

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

        self.processes = []
        """
        list: List of processes from the pipeline. This information is 
        retrieved from the .nextflow.log file in the 
        :func:`_parser_pipeline_processes` method.
        """

        self.samples = []
        """
        list: List of samples inferred from the pipeline.
        """

        self.skip_processes = ["status", "compile_status", "report",
                               "compile_reports"]
        """
        list: List of special processes that should be skipped for inspection
        purposes.
        """

        self.log_file = ".nextflow.log"

        self.pipeline_name = ""
        """
        str: Name of the pipeline, parsed from .nextflow.log
        """

        self.run_status = ""
        """
        str: Status of the pipeline. Can be either 'running', 'aborted',
        'error', 'complete'.
        """

        # CURSES ATTRIBUTES
        self.top_line = 0
        self.screen_lines = curses.LINES
        self.content_lines = 0

        self._parser_pipeline_processes()
        self._get_pipeline_status()

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
                    if process not in self.skip_processes:
                        self.processes.append(match.group(1))

                if re.match(".*Launching `.*` \[.*\] ", line):
                    match = re.match(".*Launching `.*` \[(.*)\] ", line)
                    self.pipeline_name = match.group(1)

        self.content_lines = len(self.processes)

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
        if process in self.skip_processes:
            return

        info = dict((column, fields[pos]) for column, pos in hm.items())

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

    def static_parser(self):
        """Method that parses the trace file once and updates the
        :attr:`status_info` attribute with the new entries.
        """

        # Check the timestamp of the tracefile. Only proceed with the parsing
        # if it changed from the previous time.
        timestamp = os.stat(self.trace_file)[8]
        if timestamp and timestamp == self.trace_stamp:
            return
        else:
            self.trace_stamp = timestamp

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

    def _update_process_stats(self):

        for process, vals in self.process_info.items():

            self.process_stats[process] = {}

            # Get number of samples
            self.process_stats[process]["samples"] = "{}/{}".format(
                len(vals), len(self.samples))

            # Get average time
            time_array = [self.hms(x["realtime"]) for x in vals]
            mean_time = round(sum(time_array) / len(time_array), 1)
            self.process_stats[process]["realtime"] = mean_time

            # Get cumulated time
            self.process_stats[process]["cumtime"] = sum(time_array)

    def display_overview(self):
        """Displays the default pipeline inspection overview
        """

        stay_alive = True

        screen.keypad(True)
        screen.nodelay(-1)
        curses.cbreak()
        curses.noecho()

        self.screen_lines = screen.getmaxyx()[0]

        try:
            while stay_alive:

                c = screen.getch()
                if c == curses.KEY_UP:
                    self.updown("up")
                elif c == curses.KEY_DOWN:
                    self.updown("down")
                elif c == curses.KEY_RESIZE:
                    self.screen_lines = screen.getmaxyx()[0]

                self.static_parser()
                self.flush_overview()

                sleep(0.1)
        except Exception as e:
            sys.stderr.write(e)
        finally:
            sys.stderr.write("here")
            curses.nocbreak()
            screen.keypad(0)
            curses.echo()
            curses.endwin()

    def updown(self, direction):
        """Provides curses scroll functionality.
        """

        if direction == "up" and self.top_line != 0:
            self.top_line -= 1
        elif direction == "down" and \
                screen.getmaxyx()[0] + self.top_line <= self.content_lines + 5:
            self.top_line += 1

    def flush_overview(self):
        """Displays the default overview of the pipeline execution from the
        :attr:`status_info`, :attr:`processes` and :attr:`run_status`
        attributes into stdout.
        """

        screen.erase()

        # Add static header
        screen.addstr(1, 0, "Pipeline [{}] inspection. Status: {}".format(
            self.pipeline_name, self.run_status
        ))
        screen.addstr(2, 0, "Last updated: {}".format(
            strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        headers = ["Process", "Samples", "Avg Time", "Total Time"]
        screen.addstr(4, 0, "{0: <25} "
                            "{1: ^11} "
                            "{2: ^11} "
                            "{3: ^11} ".format(*headers))

        # Get display size
        top = self.top_line
        bottom = self.screen_lines - 5 + self.top_line

        # Fetch process information
        for p, process in enumerate(self.processes[top:bottom]):

            if process not in self.process_stats:
                vals = ["-"] * 3
            else:
                ref = self.process_stats[process]
                vals = [ref["samples"], ref["realtime"], ref["cumtime"]]

            screen.addstr(
                5 + p, 0, "{0: <25} "
                          "{1: ^11} "
                          "{2: ^11} "
                          "{3: ^11} ".format(
                                process,
                                *vals
            ))

        screen.refresh()
