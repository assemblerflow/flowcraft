#!/usr/bin/env python3

import sys
import json
import logging

from os.path import dirname, abspath

logger = logging.getLogger("main.{}".format(__name__))


def write_json(report_json, version_json, trace_file, task_name,
               project_name, sample_name, pid, script_id, run_name):

    logging.info("Parsing report JSON")
    try:
        with open(report_json) as fh:
            _reports = fh.read().replace("'", '"')
            reports = json.loads(_reports)
            if "task" in reports:
                del reports["task"]
    except json.JSONDecodeError:
        logging.warning("Could not parse report JSON: {}".format(report_json))
        reports = {}

    logging.info("Parsing versions JSON")
    try:
        with open(version_json) as fh:
            _version = fh.read().replace("'", '"')
            versions = json.loads(_version)
    except json.JSONDecodeError:
        logging.warning("Could not parse versions JSON: {}".format(
            version_json))
        versions = []

    logging.info("Parsing trace file")
    with open(trace_file) as fh:
        trace = fh.readlines()

    report = {
        "pipelineId": run_name,
        "processId": pid,
        "processName": task_name,
        "projectid": run_name,
        "reportJson": reports,
        "runName": run_name,
        "scriptId": script_id,
        "versions": versions,
        "sampleName": sample_name,
        "trace": trace,
        "userId": 1,
        "username": "user",
        "workdir": dirname(abspath(report_json))
    }

    logging.info("Dumping final report JSON file")
    logging.debug("Final JSON file: {}".format(report))
    with open("{}_{}_report.json".format(task_name, sample_name), "w") \
            as report_fh:
        report_fh.write(json.dumps(report, separators=(",", ":")))


def main():

    # Fetch arguments
    args = sys.argv[1:]
    report_json = args[0]
    version_json = args[1]
    trace = args[2]
    sample_name = args[3]
    task_name = args[4]
    project_name = args[5]
    pid = args[6]
    script_id = args[7]
    run_name = args[8]
    logging.debug("Report JSON: {}".format(report_json))
    logging.debug("Version JSON: {}".format(version_json))
    logging.debug("Trace file: {}".format(trace))
    logging.debug("Sample name: {}".format(sample_name))
    logging.debug("Task name: {}".format(task_name))
    logging.debug("Project name: {}".format(project_name))
    logging.debug("Process ID: {}".format(pid))
    logging.debug("Script ID: {}".format(script_id))
    logging.debug("Run name: {}".format(run_name))

    # Write the final report JSON that compiles all information
    write_json(report_json, version_json, trace, task_name,
               project_name, sample_name, pid, script_id, run_name)


main()
