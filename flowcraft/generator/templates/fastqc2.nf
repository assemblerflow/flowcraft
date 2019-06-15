process fastqc2_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "results/fastqc2_{{ pid }}", mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}

    output:
    file "*_fastqc.{zip,html}" into {{ output_channel }}
    {% with task_name="fastqc2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    fastqc $fastq_pair
    """
}

{{ forks }}