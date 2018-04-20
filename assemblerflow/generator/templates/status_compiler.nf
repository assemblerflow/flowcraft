
/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { sample_id }
    publishDir "pipeline_status/$task_name"

    input:
    set sample_id, task_name, status, warning, fail, file(log) from {{ compile_channels }}

    output:
    file '*.status' into master_status
    file '*.warning' into master_warning
    file '*.fail' into master_fail
    file '*.log'

    """
    echo $sample_id, $task_name, \$(cat $status) > ${sample_id}_${task_name}.status
    echo $sample_id, $task_name, \$(cat $warning) > ${sample_id}_${task_name}.warning
    echo $sample_id, $task_name, \$(cat $fail) > ${sample_id}_${task_name}.fail
    echo "\$(cat .command.log)" > ${sample_id}_${task_name}.log
    """
}

process compile_status {

    publishDir 'reports/status'

    input:
    file status from master_status.collect()
    file warning from master_warning.collect()
    file fail from master_fail.collect()

    output:
    set 'master_status.csv', 'master_warning.csv', 'master_fail.csv' into mockChannel

    """
    cat $status >> master_status.csv
    cat $warning >> master_warning.csv
    cat $fail >> master_fail.csv
    """
}

