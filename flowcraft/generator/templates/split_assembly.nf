// Check for the presence of absence of the minimum contig size parameter
if (params.size{{ param_id }} == null){
    exit 1, "A minimum contig size must be provided."
}

IN_min_contig_size_{{ pid }} = Channel.value(params.size{{ param_id }})

process split_assembly_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/assembly/split_assembly_{{ pid }}/${sample_id}/"

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    val min_contig_size from IN_min_contig_size_{{ pid }}

    output:
    file '*split.fasta' into splitCh_{{ pid }} optional true
    {% with task_name="split_assembly" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "split_fasta.py"


}
{{ output_channel }} = Channel.create()

splitCh_{{ pid }}.flatMap().map{ it -> [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it] }
    .into({{ output_channel }})

{{ forks }}