
/** Reports
Compiles the reports from every process
*/
process report {

    tag { sample_id }

    input:
    set sample_id,
            task_name,
            pid,
            report_json,
            version_json,
            trace from {{ compile_channels }}

    output:
    file "*" optional true into master_report

    """
    prepare_reports.py $report_json $version_json $trace $sample_id $task_name 1 $pid
    """

}


process compile_reports {

    publishDir "pipeline_report/"

    input:
    file report from master_report.collect()

    output:
    file "pipeline_report.json"

    """
    #!/usr/bin/env python3
    import sys
    import json

    reports = '${report}'.split()

    storage = []
    for r in reports:
        with open(r) as fh:
            rjson = json.load(fh)
            storage.append(rjson)
            print("{}: {}".format(rjson["processName"], sys.getsizeof(json.dumps(rjson))))

    with open("pipeline_report.json", "w") as rep_fh:
       rep_fh.write(json.dumps({"data": {"results": storage}},
                    separators=(",", ":")))
    """

}

