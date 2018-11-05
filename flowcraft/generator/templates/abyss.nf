process abyss_{{ pid }} {
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/abyss_{{ pid }}/', pattern: '*-scaffolds.fa'
    publishDir 'results/assembly/abyss_{{ pid }}/', pattern: '*-scaffolds.gfa'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val k from Channel.value(params.abyssKmer{{ param_id }})

    output:
    set sample_id, file('*-scaffolds.fa') into {{ output_channel }}
    file "*-scaffolds.gfa" into gfa1_{{ pid }}
    {% with task_name="abyss" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    "abyss-pe name=${sample_id} graph=gfa k=${k} v=-v in=\"${fastq_pair[0]} ${fastq_pair[1]}\""
}

{{ forks }}
