process metaspades_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }
    publishDir 'results/assembly/spades_{{ pid }}/', pattern: '*_metaspades*.fasta', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), max_len from {{ input_channel }}.join(SIDE_max_len_{{ pid }})
    val kmers from IN_metaspades_kmers

    output:
    set fastq_id, file('*_metaspades*.fasta') into {{ output_channel }}
    {% with task_name="metaspades" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "metaspades.py"

}

{{ forks }}