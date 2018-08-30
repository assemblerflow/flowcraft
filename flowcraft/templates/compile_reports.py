#!/usr/bin/python3
import sys
import json
import logging

REPORTS = "${report}".split()
FORKS = "${forks}"
DAG = "${dag}"


def main(reports, forks, dag):

    metadata = {
        "nfMetadata": {
            "scriptId": "${workflow.scriptId}",
            "scriptName": "${workflow.scriptId}",
            "profile": "${workflow.profile}",
            "container": "${workflow.container}",
            "containerEngine": "${workflow.containerEngine}",
            "commandLine": "${workflow.commandLine}",
            "runName": "${workflow.runName}",
            "sessionId": "${workflow.sessionId}",
            "projectDir": "${workflow.projectDir}",
            "launchDir": "${workflow.launchDir}",
            "start_time": "${workflow.start}"
        }
    }

    # Add nextflow metadata
    storage = []

    # Add forks dictionary
    try:
        with open(forks) as fh:
            forks = json.load(fh)
            metadata["nfMetadata"]["forks"] = forks
    except json.JSONDecodeError:
        logging.warning("Could not parse versions JSON: {}".format(
            dag))

    # Add tree DAG in JSON format
    try:
        with open(dag) as fh:
            dag = json.load(fh)
            metadata["nfMetadata"]["dag"] = dag
    except json.JSONDecodeError:
        logging.warning("Could not parse versions JSON: {}".format(
            dag))

    storage.append(metadata)
    # Write metadata information to dotfile. This dotfile is then sent to the
    # ReportHTTP, when available in the afterScript process directive.
    with open(".metadata.json", "w") as fh:
        fh.write(json.dumps(metadata, separators=(",", ":")))

    for r in reports:
        with open(r) as fh:
            rjson = json.load(fh)
            storage.append(rjson)
            print("{}: {}".format(rjson["processName"],
                                  sys.getsizeof(json.dumps(rjson))))

    with open("pipeline_report.json", "w") as rep_fh:
        rep_fh.write(json.dumps({"data": {"results": storage}},
                                separators=(",", ":")))


if __name__ == "__main__":
    main(REPORTS, FORKS, DAG)
