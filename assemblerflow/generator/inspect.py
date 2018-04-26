import re
import os

from os.path import join, abspath


class NextflowInspector:

    def __init__(self, trace_file):

        self.trace_file = trace_file
        """
        str: Path to nextflow trace file.
        """

        self.stored_ids = []
        """
        list: Stores the task_ids that have already been parsed
        """

        self.status_info = {}
        """
        dict: Main object that stores the status information for each process
        name in the trace file.
        """

        self.processes = []
        """
        list: List of processes from the pipeline. This information is 
        retrieved from the .nextflow.log file in the 
        :func:`_parser_pipeline_processes` method.
        """

        self.skip_processes = ["status", "compile_status", "report",
                               "compile_reports"]
        """
        list: List of special processes that should be skipped for inspection
        purposes.
        """

        self.log_file = ".nextflow.log"

        self.run_status = ""
        """
        str: Status of the pipeline. Can be either 'running', 'aborted',
        'error', 'complete'.
        """

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

        self.status_info[process] = \
            dict((column, fields[pos]) for column, pos in hm.items())

        # If the task hash code is provided, expand it to the work directory
        # and add a new entry
        if self.status_info[process]["hash"]:
            hs = self.status_info[process]["hash"]
            self.status_info[process]["work_dir"] = self._expand_path(hs)

    def static_parser(self):
        """Method that parses the trace file once and updates the
        :attr:`status_info` attribute with the new entries.
        """

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

