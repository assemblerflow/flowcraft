// Check sliding window parameter
if ( params.trimSlidingWindow{{ param_id }}.toString().split(":").size() != 2 ){
    exit 1, "'trimSlidingWindow{{ param_id}}' parameter must contain two values separated by a ':'. Provided value: '${params.trimSlidingWindow{{ param_id}}}'"
}
if ( !params.trimLeading{{ param_id}}.toString().isNumber() ){
    exit 1, "'trimLeading{{ param_id}}' parameter must be a number. Provide value: '${params.trimLeading_{{pid}}}'"
}
if ( !params.trimTrailing{{ param_id}}.toString().isNumber() ){
    exit 1, "'trimTrailing{{ param_id}}' parameter must be a number. Provide value: '${params.trimTrailing{{ param_id}}}'"
}
if ( !params.trimMinLength{{ param_id}}.toString().isNumber() ){
    exit 1, "'trimMinLength{{ param_id}}' parameter must be a number. Provide value: '${params.trimMinLength{{ param_id}}}'"
}

IN_trimmomatic_opts_{{ pid }} = Channel.value([params.trimSlidingWindow{{ param_id}},params.trimLeading{{ param_id}},params.trimTrailing{{ param_id}},params.trimMinLength{{ param_id}}])
IN_adapters_{{ pid }} = Channel.value(params.adapters{{ param_id}})

clear = params.clearAtCheckpoint ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process trimmomatic_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    publishDir "results/trimmomatic_{{ pid }}", pattern: "*.gz"

    tag { sample_id }

    input:
    set sample_id, file(fastq_pair), phred from {{ input_channel }}.join(SIDE_phred_{{ pid }})
    val trim_range from Channel.value("None")
    val opts from IN_trimmomatic_opts_{{ pid }}
    val ad from IN_adapters_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, "${sample_id}_*trim.fastq.gz" into {{ output_channel }}
    file 'trimmomatic_report.csv'
    {% with task_name="trimmomatic" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "trimmomatic.py"

}

{{ forks }}

