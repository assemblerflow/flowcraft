
/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { fastq_id }
    publishDir "pipeline_status/$task_name"

    input:
    set fastq_id, task_name, status, warning, fail, file(log) from {{ status_channels }}

    output:
    file '*.status' into master_status
    file '*.warning' into master_warning
    file '*.fail' into master_fail
    file '*.log'

    """
    echo $fastq_id, $task_name, \$(cat $status) > ${fastq_id}.status
    echo $fastq_id, $task_name, \$(cat $warning) > ${fastq_id}.warning
    echo $fastq_id, $task_name, \$(cat $fail) > ${fastq_id}.fail
    echo "\$(cat .command.log)" > ${fastq_id}.log
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

