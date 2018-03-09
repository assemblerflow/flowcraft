
process process_spades {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    // This process can only use a single CPU
    cpus 1
    publishDir "reports/assembly/spades_filter", pattern: '*.report.csv', mode: 'copy'

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    val opts from IN_process_spades_opts
    val gsize from IN_genome_size

    output:
    set fastq_id, file('*.assembly.fasta') optional true into {{ output_channel }}
    file '*.report.csv' optional true
    {% with task_name="process_spades" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    when:
    params.stopAt != "process_spades"

    script:
    template "process_spades.py"

}

{{ forks }}

