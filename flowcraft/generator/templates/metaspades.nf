if ( params.metaspadesKmers{{ param_id }}.toString().split(" ").size() <= 1 ){
    if (params.metaspadesKmers{{ param_id }}.toString() != 'auto'){
        exit 1, "'metaspadesKmers{{ param_id }}' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.metaspadesKmers{{ param_id }}}"
    }
}
IN_metaspades_kmers_{{pid}} = Channel.value(params.metaspadesKmers{{ param_id }})


process metaspades_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/metaspades_{{ pid }}/', pattern: '*_metaspades*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from {{ input_channel }}.join(SIDE_max_len_{{ pid }})
    val kmers from IN_metaspades_kmers_{{pid}}

    output:
    set sample_id, file('*_metaspades*.fasta') into {{ output_channel }}
    {% with task_name="metaspades" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "metaspades.py"

}

{{ forks }}