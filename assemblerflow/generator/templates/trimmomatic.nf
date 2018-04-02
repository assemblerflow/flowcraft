
process trimmomatic_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }

    input:
    set fastq_id, file(fastq_pair), phred from {{ input_channel }}.join(SIDE_phred_{{ pid }})
    val trim_range from Channel.value("None")
    val opts from IN_trimmomatic_opts
    val ad from IN_adapters

    output:
    set fastq_id, "${fastq_id}_*P*" optional true into {{ output_channel }}
    file 'trimmomatic_report.csv'
    {% with task_name="trimmomatic" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    when:
    params.stopAt != "trimmomatic"

    script:
    template "trimmomatic.py"

}

{{ forks }}

