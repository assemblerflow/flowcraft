
process spades_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/spades_{{ pid }}/', pattern: '*_spades*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from {{ input_channel }}.join(SIDE_max_len_{{ pid }})
    val opts from IN_spades_opts
    val kmers from IN_spades_kmers

    output:
    set sample_id, file('*_spades*.fasta') into {{ output_channel }}
    {% with task_name="spades" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "spades.py"

}

{{ forks }}

