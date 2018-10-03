process metaphlan_fa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/metaphlan/", pattern: "*.txt"

    input:
    set sample_id, file(fasta) from {{ input_channel }}

    output:
     set sample_id, file("${sample_id}_profiled_metagenome.txt") into {{ output_channel }}
    {% with task_name="metaphlan_fa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    metaphlan2.py ${fasta} --nproc $task.cpus --input_type fasta > ${sample_id}_profiled_metagenome.txt

    """
}

{{ forks }}