
/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { fastq_id }

    input:
    set fastq_id, task_name, status, warning, fail from {{ status_channels }}

    output:
    file 'status_*' into master_status
    file 'warning_*' into master_warning
    file 'fail_*' into master_fail

    """
    echo $fastq_id, $task_name, \$(cat $status) > status_${fastq_id}_${task_name}
    echo $fastq_id, $task_name, \$(cat $warning) > warning_${fastq_id}_${task_name}
    echo $fastq_id, $task_name, \$(cat $fail) > fail_${fastq_id}_${task_name}
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

