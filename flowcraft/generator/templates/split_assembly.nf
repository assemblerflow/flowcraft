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
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})
    val min_contig_size from IN_min_contig_size_{{ pid }}

    output:
    set file('*.fasta'), file(fastq) into splitCh_{{ pid }}
    {% with task_name="split_assembly" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "split_fasta.py"


}

{{ output_channel }} = Channel.create()

splitCh_{{ pid }}.map{ it -> [it[0].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it[0], it[1]]}.into({{ output_channel }})

{{ forks }}