#!/usr/bin/env python3

import sys
import json

from os.path import dirname, abspath


def write_json(report_json, task_name, sample_name, pid):

    with open(report_json) as fh:
        res = json.load(fh)

    res["task"] = task_name
    del res["task"]

    report = {
        "reportJson": res,
        "processId": pid,
        "pipelineId": 1,
        "projectid": 1,
        "userId": 1,
        "username": "user",
        "processName": task_name,
        "workdir": dirname(abspath(report_json))
    }

    with open("{}_{}_report.json".format(task_name, sample_name), "w") \
            as report_fh:
        report_fh.write(json.dumps(report, separators=(",", ":")))


def main():

    args = sys.argv[1:]
    report_json = args[0]
    sample_name = args[1]
    task_name = args[2]
    project_name = args[3]
    pid = args[4]

    print(report_json)
    print(sample_name)
    print(task_name)
    print(project_name)
    print(pid)

    try:
        write_json(report_json, task_name, sample_name, pid)
    except json.decoder.JSONDecodeError:
        print("Could not parse JSON output from {}, sample name {} and "
              "pid {}".format(report_json, sample_name, pid))


main()
