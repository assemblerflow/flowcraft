
process trimmomatic_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    publishDir "results/trimmomatic_{{ pid }}", pattern: "*.gz"

    tag { sample_id }

    input:
    set sample_id, file(fastq_pair), phred from {{ input_channel }}.join(SIDE_phred_{{ pid }})
    val trim_range from Channel.value("None")
    val opts from IN_trimmomatic_opts
    val ad from IN_adapters

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

