if ( !params.skesaMinKmerCoverage{{ param_id }}.toString().isNumber() ){ 
    exit 1, "'skesaMinKmerCoverage{{ param_id }}' parameter must be a number. Provided value: ${params.skesaMinKmerCoverage{{ param_id }}}"
}
if ( !params.skesaMinContigLen{{ param_id }}.toString().isNumber() ){ 
    exit 1, "'skesaMinContigLen{{ param_id }}' parameter must be a number. Provided value: ${params.skesaMinContigLen{{ param_id }}}"
}
if ( !params.skesaMaxContigs{{ param_id }}.toString().isNumber() ){ 
    exit 1, "'skesaMaxContigs{{ param_id }}' parameter must be a number. Provided value: ${params.skesaMaxContigs{{ param_id }}}"
}

IN_process_skesa_opts_{{ pid }} = Channel.value([params.skesaMinContigLen{{ param_id }},params.skesaMinKmerCoverage{{ param_id }},params.skesaMaxContigs{{ param_id }}])
IN_genome_size_{{ pid }} = Channel.value(params.genomeSize{{ param_id }})

process process_skesa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir "reports/assembly/skesa_filter_{{ pid }}", pattern: '*.report.csv', mode: 'copy'

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    val opts from IN_process_skesa_opts_{{ pid }}
    val gsize from IN_genome_size_{{ pid }}
    val assembler from Channel.value("skesa")

    output:
    set sample_id, file('*.fasta') into {{ output_channel }}
    file '*.report.csv' optional true
    {% with task_name="process_skesa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_assembly.py"

}

{{ forks }}
