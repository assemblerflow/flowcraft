
/** Reports
Compiles the reports from every process
*/
process report {

    input:
    set fastq_id, task_name, pid, report_json from {{ compile_channels }}

    output:
    file "*.json" into master_report

    """
    prepare_reports.py $report_json $fastq_id $task_name 1 $pid
    """

}


process compile_reports {

    publishDir "reports/json_reports"

    input:
    file report from master_report.collect()

    output:
    file "*.json"

    """

    """

}

