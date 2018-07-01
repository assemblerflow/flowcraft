if ( params.megahitKmers{{ param_id }}.toString().split(" ").size() <= 1 ){
    if (params.megahitKmers{{ param_id }}.toString() != 'auto'){
        exit 1, "'megahitKmers{{ param_id }}' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.megahitKmers{{ param_id }}}"
    }
}
IN_megahit_kmers_{{ pid }} = Channel.value(params.megahitKmers{{ param_id }})


process megahit_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/megahit_{{ pid }}/', pattern: '*_megahit*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from {{ input_channel }}.join(SIDE_max_len_{{ pid }})
    val kmers from IN_megahit_kmers_{{ pid }}

    output:
    set sample_id, file('*megahit*.fasta') into {{ output_channel }}
    {% with task_name="megahit" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "megahit.py"

}

{{ forks }}



