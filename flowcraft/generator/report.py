import os
import re
import sys
import json
import signal
import socket
import hashlib
import logging
import requests

from os.path import join, abspath
from time import sleep
from pympler.asizeof import asizeof

try:
    import generator.error_handling as eh
    from generator.process_details import colored_print
except ImportError:
    import flowcraft.generator.error_handling as eh
    from flowcraft.generator.process_details import colored_print

logger = logging.getLogger("main.{}".format(__name__))


def signal_handler():
    """This function is bound to the SIGINT signal (like ctrl+c) to graciously
    exit the program and reset the curses options.
    """

    print("Exiting flowcraft report brodcast... Bye")
    sys.exit(0)


class FlowcraftReport:

    def __init__(self, report_file, trace_file=None, log_file=None,
                 watch=False, ip_addr=None):

        self.report_file = report_file
        """
        str: Path to Report JSON file.
        """

        if not ip_addr:
            self.app_address = "http://192.92.149.169:80/"
        else:
            self.app_address = ip_addr
            """
            str: Address of flowcraft web app
            """

        self.broadcast_address = "{}reports/broadcast/api/reports".format(
            self.app_address)

        self.refresh_rate = 1

        self.send = True
        """
        boolean: This attribute is used when the report mode is used with the
        --watch option. It will be set to False after sending a request, and 
        set to True when there is a change in the pipeline reports.
        """

        self.watch = watch
        """
        boolean: When False, the reports mode will try to open the provided
        report JSON file and send it to the flowcraft service. When True, 
        it will try to open the nextflow trace file instead and continuously 
        compile the report JSON files from the `report` processes as they 
        are created. 
        """

        self.log_file = log_file
        """
        str: Path to .nextflow.log file.
        """

        self.log_sizestamp = None
        """
        str: Stores the sizestamp of the last modification of the trace file.
        This is used to parse the file only when it has changed.
        """

        self.status_info = None
        """
        str: Status of the pipeline execution. Used in the watch report mode
        and varies between 'running', 'aborted', 'complete'.
        """

        self.trace_file = trace_file
        """
        str: Path to nextflow trace file.
        """

        self.trace_sizestamp = None
        """
        str: Stores the sizestamp of the last modification of the trace file.
        This is used to parse the file only when it has changed.
        """

        self.trace_retry = 0
        """
        int: Each time the log file is not found, this counter is 
        increased. Only when it matches the :attr:`MAX_RETRIES` attribute
        does it raises a FileNotFoundError.
        """

        self.stored_ids = []
        """
        list: Stores the task_ids that have already been parsed. It is used
        to skip them when parsing the trace files multiple times.
        """

        self.report_queue = []
        """
        list: Stores the paths of the report JSON files that are on queue to
        be sent to the flowcraft service. This list will be emptied when these
        JSONs are sent.
        """

        # Checks if report file is available
        self._check_required_files()

        signal.signal(signal.SIGINT, lambda *x: signal_handler())

    def _check_required_files(self):

        if not os.path.exists(self.report_file) and not self.watch:
            raise eh.ReportError("The provided report JSON file could not be"
                                 " opened: {}".format(self.report_file))

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

    def _get_report_id(self):
        """Returns a hash of the reports JSON file
        """

        if self.watch:

            with open(self.log_file) as fh:
                header = fh.readline()

            pipeline_path = re.match(
                ".*nextflow run ([^\s]+).*", header).group(1)

            # Get hash from the entire pipeline file
            pipeline_hash = hashlib.md5()
            with open(pipeline_path, "rb") as fh:
                for chunk in iter(lambda: fh.read(4096), b""):
                    pipeline_hash.update(chunk)
            # Get hash from the current working dir and hostname
            workdir = os.getcwd().encode("utf8")
            hostname = socket.gethostname().encode("utf8")
            dir_hash = hashlib.md5(workdir + hostname)

            return pipeline_hash.hexdigest() + dir_hash.hexdigest()

        else:
            with open(self.report_file) as fh:
                report_json = json.loads(fh.read())

            metadata = report_json["data"]["results"][0]["nfMetadata"]

            try:
                report_id = metadata["scriptId"] + metadata["sessionId"]
            except KeyError:
                raise eh.ReportError("Incomplete or corrupt report JSON file "
                                     "missing the 'scriptId' and/or 'sessionId' "
                                     "metadata information")

            return report_id

    def _update_pipeline_status(self):
        """
        Parses the .nextflow.log file for signatures of pipeline status and sets
        the :attr:`status_info` attribute.
        """

        prev_status = self.status_info

        with open(self.log_file) as fh:

            for line in fh:

                if "Session aborted" in line:
                    self.status_info = "aborted"
                    self.send = True if prev_status != self.status_info \
                        else self.send
                    return

                if "Execution complete -- Goodbye" in line:
                    self.status_info = "complete"
                    self.send = True if prev_status != self.status_info \
                        else self.send
                    return

            self.status_info = "running"
            self.send = True if prev_status != self.status_info \
                else self.send

    def update_trace_watch(self):
        """Parses the nextflow trace file and retrieves the path of report JSON
        files that have not been sent to the service yet.
        """

        # Check the size stamp of the tracefile. Only proceed with the parsing
        # if it changed from the previous size.
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

                if fields[hm["process"]] == "report":
                    self.report_queue.append(
                        self._expand_path(fields[hm["hash"]])
                    )
                    self.send = True

                # Add the processed trace line to the stored ids. It will be
                # skipped in future parsers
                self.stored_ids.append(fields[hm["task_id"]])

    def update_log_watch(self):
        """Parses nextflow log file and updates the run status
        """

        # Check the size stamp of the tracefile. Only proceed with the parsing
        # if it changed from the previous size.
        size_stamp = os.path.getsize(self.log_file)
        self.trace_retry = 0
        if size_stamp and size_stamp == self.log_sizestamp:
            return
        else:
            logger.debug("Updating log size stamp to: {}".format(size_stamp))
            self.log_sizestamp = size_stamp

        self._update_pipeline_status()

    def _send_live_report(self, report_id):
        """Sends a PUT request with the report JSON files currently in the
        report_queue attribute.

        Parameters
        ----------
        report_id : str
            Hash of the report JSON as retrieved from :func:`~_get_report_hash`
        """

        # Determines the maximum number of reports sent at the same time in
        # the same payload
        buffer_size = 100
        logger.debug("Report buffer size set to: {}".format(buffer_size))

        for i in range(0, len(self.report_queue), buffer_size):

            # Reset the report compilation batch
            reports_compilation = []

            # Iterate over report JSON batches determined by buffer_size
            for report in self.report_queue[i: i + buffer_size]:
                try:
                    report_file = [x for x in os.listdir(report)
                                   if x.endswith(".json")][0]
                except IndexError:
                    continue
                with open(join(report, report_file)) as fh:
                    reports_compilation.append(json.loads(fh.read()))

            logger.debug("Payload sent with size: {}".format(
                asizeof(json.dumps(reports_compilation))
            ))
            logger.debug("status: {}".format(self.status_info))

            try:
                requests.put(
                    self.broadcast_address,
                    json={"run_id": report_id,
                          "report_json": reports_compilation,
                          "status": self.status_info}
                )
            except requests.exceptions.ConnectionError:
                logger.error(colored_print(
                    "ERROR: Could not establish connection with server. The server"
                    " may be down or there is a problem with your internet "
                    "connection.", "red_bold"))
                sys.exit(1)

        # When there is no change in the report queue, but there is a change
        # in the run status of the pipeline
        if not self.report_queue:

            logger.debug("status: {}".format(self.status_info))

            try:
                requests.put(
                    self.broadcast_address,
                    json={"run_id": report_id,
                          "report_json": [],
                          "status": self.status_info}
                )
            except requests.exceptions.ConnectionError:
                logger.error(colored_print(
                    "ERROR: Could not establish connection with server. The"
                    " server may be down or there is a problem with your "
                    "internet connection.", "red_bold"))
                sys.exit(1)

        # Reset the report queue after sending the request
        self.report_queue = []

    def _init_live_reports(self, report_id):
        """Sends a POST request to initialize the live reports

        Parameters
        ----------
        report_id : str
            Hash of the report JSON as retrieved from :func:`~_get_report_hash`
        """

        logger.debug("Sending initial POST request to {} to start report live"
                     " update".format(self.broadcast_address))

        try:
            with open(".metadata.json") as fh:
                metadata = [json.load(fh)]
        except:
            metadata = []

        start_json = {
            "data": {"results": metadata}
        }

        try:
            requests.post(
                self.broadcast_address,
                json={"run_id": report_id, "report_json": start_json,
                      "status": self.status_info}
            )
        except requests.exceptions.ConnectionError:
            logger.error(colored_print(
                "ERROR: Could not establish connection with server. The server"
                " may be down or there is a problem with your internet "
                "connection.", "red_bold"))
            sys.exit(1)

    def _close_connection(self, report_id):
        """Sends a delete request for the report JSON hash

        Parameters
        ----------
        report_id : str
            Hash of the report JSON as retrieved from :func:`~_get_report_hash`
        """

        logger.debug(
            "Closing connection and sending DELETE request to {}".format(
                self.broadcast_address))

        try:
            r = requests.delete(self.broadcast_address,
                                json={"run_id": report_id})
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

    def _send_report(self, report_id):

        with open(self.report_file) as fh:
            report_json = json.loads(fh.read())

        logger.debug("Unique payload sent with size: {}".format(
            asizeof(json.dumps(report_json))
        ))

        try:
            requests.post(
                self.broadcast_address,
                json={"run_id": report_id, "report_json": report_json}
            )
        except requests.exceptions.ConnectionError:
            logger.error(colored_print(
                "ERROR: Could not establish connection with server. The server"
                " may be down or there is a problem with your internet "
                "connection.", "red_bold"))
            sys.exit(1)

    def _print_msg(self, run_id):

        report_address = "{}reports/broadcast/{}".format(self.app_address,
                                                         run_id)
        logger.info(colored_print(
            "The pipeline reports are available in the following link:",
            "green_bold"))
        logger.info("{}".format(report_address))

    def broadcast_report(self):

        logger.info(colored_print("Preparing to broacast reports...",
                                  "green_bold"))

        report_hash = self._get_report_id()

        # When in watch mode,
        if self.watch:
            logger.info(colored_print("\tFetching pipeline run status",
                                      "green_bold"))
            self._update_pipeline_status()
            logger.info(colored_print(
                "\tSending initial request to test service", "green_bold"))
            self._init_live_reports(report_hash)
            logger.info(colored_print("\tInitial parsing of trace file",
                                      "green_bold"))
            self.update_trace_watch()

            self._print_msg(report_hash)

        logger.debug("Establishing connection...")

        stay_alive = True
        _broadcast_sent = False
        try:
            while stay_alive:

                # When not in watch mode, send the report JSON once
                if not _broadcast_sent and not self.watch:
                    self._send_report(report_hash)
                    self._print_msg(report_hash)
                    _broadcast_sent = True

                # When in watch mode, continuously monitor the trace file for
                # updates
                if self.watch:
                    self.update_trace_watch()
                    self.update_log_watch()
                    # When new report JSON files are available, send then
                    # via a PUT request
                    if self.send:
                        self._send_live_report(report_hash)
                        self.send = False

                sleep(self.refresh_rate)

        except FileNotFoundError as e:
            print(e)
            logger.error(colored_print(
                "ERROR: Report JSON file is not reachable!", "red_bold"))
        except Exception as e:
            logger.exception("ERROR: " + e)
        finally:
            logger.info("Closing connection")
            self._close_connection(report_hash)
