#!/usr/bin/env python

import sys
import json


def write_json(report_json, task_name, project_name, sample_name, pid):

    with open(report_json) as fh:
        res = json.load(fh)

    report = {
        "report_json": res,
        "task": task_name,
        "project_id": project_name,
        "sample_name": sample_name,
        "pipeline_id": sample_name,
        "process_id": pid
    }

    with open("{}_{}_report.json".format(task_name, pid), "w") as report_fh:
        report_fh.write(json.dumps(report, separators=(",", ":")))


def main():

    args = sys.argv[1:]
    report_json = args[0]
    sample_name = args[1]
    task_name = args[2]
    project_name = args[3]
    pid = args[4]

    write_json(report_json, task_name, project_name, sample_name, pid)

    print(report_json)
    print(task_name)
    print(project_name)


main()
