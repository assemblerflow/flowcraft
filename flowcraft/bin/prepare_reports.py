#!/usr/bin/env python3

import sys
import json


def write_json(report_json, task_name, project_name, sample_name, pid):

    with open(report_json) as fh:
        res = json.load(fh)

    res["task"] = task_name

    report = {
        "report_json": res,
        "project_id": project_name,
        "sample_name": sample_name,
        "pipeline_id": sample_name,
        "process_id": pid,
        "user_id": 1,
        "username": "user"
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
        write_json(report_json, task_name, project_name, sample_name, pid)
    except json.decoder.JSONDecodeError:
        print("Could not parse JSON output from {}, sample name {} and "
              "pid {}".format(report_json, sample_name, pid))


main()
