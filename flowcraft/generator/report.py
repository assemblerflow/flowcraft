import os
import sys
import json
import signal
import logging
import requests

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

    def __init__(self, report_file, ip_addr=None):

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

        # Checks if report file is available
        self._check_required_files()

        signal.signal(signal.SIGINT, lambda *x: signal_handler())

    def _check_required_files(self):
        if not os.path.exists(self.report_file):
            raise eh.ReportError("The provided report JSON file could not be"
                                 " opened: {}".format(self.report_file))

    def _get_report_id(self):
        """Returns a hash of the reports JSON file
        """

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

    def _close_connection(self, report_id):
        """Sends a delete request for the report JSON hash

        Parameters
        ----------
        report_id : str
            Hash of the report JSON as retrieved from :func:`~_get_report_hash`
        """

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

        logger.debug("Payload sent with size: {}".format(
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

        logger.debug("Establishing connection...")

        stay_alive = True
        _broadcast_sent = False
        try:
            while stay_alive:

                if not _broadcast_sent:
                    self._send_report(report_hash)
                    self._print_msg(report_hash)
                    _broadcast_sent = True

                sleep(self.refresh_rate)

        except FileNotFoundError:
            logger.error(colored_print(
                "ERROR: Report JSON file is not reachable!", "red_bold"))
        except Exception as e:
            logger.exception("ERROR: " + e)
        finally:
            logger.info("Closing connection")
            self._close_connection(report_hash)
