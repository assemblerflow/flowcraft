
clear = params.clearAtCheckpoint ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process skesa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/skesa_{{ pid }}', pattern: '*skesa*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, file('*.fasta') into {{ output_channel }}
    {% with task_name="skesa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "skesa.py"

}

{{ forks }}