#!/usr/bin/python3
import os
import sys
import json
import zipfile
import logging

REPORTS = "${report}".split()
FORKS = "${forks}"
DAG = "${dag}"
MAIN_JS = "${js}"


html_template = """
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <link href="https://fonts.googleapis.com/icon?family=Material+Icons" rel="stylesheet">
  <title>FlowCraft App</title>
</head>
<body style="background-color: #f2f2f2">
    <div id="app"><!-- React --></div>
</body>
<script> const _fileReportData = {} </script>
<script src="./src/main.js"></script>
</html>
"""


def main(reports, forks, dag, main_js):

    metadata = {
        "nfMetadata": {
            "scriptId": "${workflow.scriptId}",
            "scriptName": "${workflow.scriptName}",
            "profile": "${workflow.profile}",
            "container": "${workflow.container}",
            "containerEngine": "${workflow.containerEngine}",
            "commandLine": "${workflow.commandLine}",
            "runName": "${workflow.runName}",
            "sessionId": "${workflow.sessionId}",
            "projectDir": "${workflow.projectDir}",
            "launchDir": "${workflow.launchDir}",
            "startTime": "${workflow.start}"
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

    with open("pipeline_report.html", "w") as html_fh:
        html_fh.write(html_template.format(
            json.dumps({"data": {"results": storage}}, separators=(",", ":"))))

    with zipfile.ZipFile(MAIN_JS) as zf:
        os.mkdir("src")
        zf.extractall("./src")

    with open("pipeline_report.json", "w") as rep_fh:
        rep_fh.write(json.dumps({"data": {"results": storage}},
                                separators=(",", ":")))


if __name__ == "__main__":
    main(REPORTS, FORKS, DAG, MAIN_JS)
