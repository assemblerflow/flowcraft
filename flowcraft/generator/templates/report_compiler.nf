
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

    publishDir "pipeline_report/", mode: "copy"

    input:
    file report from master_report.collect()
    file forks from Channel.fromPath(".forkTree.json")
    file dag from Channel.fromPath(".treeDag.json")

    output:
    file "pipeline_report.json"

    """
    #!/usr/bin/env python3
    import sys
    import json

    reports = "${report}".split()
    forks = "${forks}"
    dag = "${dag}"

    # Add nextflow metadata
    storage = [{
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
    }]

    # Add forks dictionary
    try:
        with open(forks) as fh:
            forks = json.load(fh)
            storage[0]["nfMetadata"]["forks"] =  forks
    except json.JSONDecodeError:
        logging.warning("Could not parse versions JSON: {}".format(
            dag))

    # Add tree DAG in JSON format
    try:
        with open(dag) as fh:
            dag = json.load(fh)
            storage[0]["nfMetadata"]["dag"] =  dag
    except json.JSONDecodeError:
        logging.warning("Could not parse versions JSON: {}".format(
            dag))

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

