if ( !params.spadesMinKmerCoverage{{ param_id }}.toString().isNumber()){
    exit 1, "'spadesMinKmerCoverage' parameter must be a number. Provided value: ${params.spadesMinKmerCoverage{{ param_id }}}"
}
if ( !params.spadesMinContigLen{{ param_id }}.toString().isNumber() ){
    exit 1, "'spadesMinContigLen' parameter must be a number. Provided value: ${params.spadesMinContigLen{{ param_id }}}"
}
if ( !params.spadesMaxContigs{{ param_id }}.toString().isNumber() ){
    exit 1, "'spadesMaxContigs' parameter must be a number. Provided value: ${params.spadesMaxContigs{{ param_id }}}"
}

IN_process_spades_opts_{{ pid }} = Channel.value([params.spadesMinContigLen{{ param_id }}, params.spadesMinKmerCoverage{{ param_id }}, params.spadesMaxContigs{{ param_id }}])
IN_genome_size_{{ pid }} = Channel.value(params.genomeSize{{ param_id }})

process process_spades_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir "reports/assembly/spades_filter_{{ pid }}", pattern: '*.report.csv', mode: 'copy'

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    val opts from IN_process_spades_opts_{{ pid }}
    val gsize from IN_genome_size_{{ pid }}
    val assembler from Channel.value("spades")

    output:
    set sample_id, file('*.fasta') into {{ output_channel }}
    file '*.report.csv' optional true
    {% with task_name="process_spades" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_assembly.py"

}

{{ forks }}

