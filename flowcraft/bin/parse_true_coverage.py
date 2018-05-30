#!/usr/bin/env python

import sys
import json


def parse_true_coverage(report_json, fail_json=None):

    with open(report_json) as fh:
        res = json.load(fh)
        print("Report JSON: {}".format(res))

    with open(".report.json", "w") as report_fh:

        json_dic = {
            "tableRow": [
                {"header": "True Coverage",
                 "value": res["mean_sample_coverage"],
                 "table": "assembly",
                 "columnBar": True},
            ]
        }

        if fail_json:
            with open(fail_json) as fail_fh:
                fail = json.load(fail_fh)
                print("Fail JSON: {}".format(fail))

            json_dic["fail"] = {
                "process": "true_coverage",
                "value": []
            }

            for v in fail.values():
                json_dic["fail"]["value"].append(v)

        report_fh.write(json.dumps(json_dic, separators=(",", ":")))


def main():

    args = sys.argv[1:]
    report_json = args[0]
    try:
        fail_json = args[1]
    except IndexError:
        fail_json = None

    print("Parsing report {} and fail {}".format(report_json, fail_json))

    parse_true_coverage(report_json, fail_json)


main()
