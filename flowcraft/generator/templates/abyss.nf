
IN_abyss_kmer_{{ pid }} = Channel.value(params.abyssKmer{{ param_id }})

process abyss_{{ pid }} {
    
    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/abyss_{{ pid }}/', pattern: '*-scaffolds.fa', mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val kmer from IN_abyss_kmer_{{ pid }}

    output:
    set sample_id, file('*.fa') into {{ output_channel }}
    {% with task_name="abyss" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    "abyss-pe in=\"${fastq_pair[0]} ${fastq_pair[1]}\" k=${kmer} name=${sample_id}"
}

{{ forks }}
