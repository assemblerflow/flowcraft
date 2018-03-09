
process spades {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    publishDir 'results/assembly/spades/', pattern: '*_spades.assembly.fasta', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), max_len from {{ input_channel }}.join(SIDE_max_len_{{ pid }})
    val opts from IN_spades_opts
    val kmers from IN_spades_kmers

    output:
    set fastq_id, file('*_spades.assembly.fasta') optional true into {{ output_channel }}
    {% with task_name="spades" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    when:
    params.stopAt != "spades"

    script:
    template "spades.py"

}

{{ forks }}

