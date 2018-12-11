// True when a dengue_typing secondary channel is connected
has_ref_{{pid}} = binding.hasVariable('_ref_seqTyping_{{ pid }}')

if ( has_ref_{{pid}} ){
    {{ input_channel }}.map{ it[1] }.collect().mix(_ref_seqTyping_{{pid}}.unique{it.name}).set{mafft_input}
} else {
    {{ input_channel }}.map{ it[1] }.collect().set{mafft_input}
}

//{{ input_channel }}.map{ it[1] }.mix(_ref_seqTyping_{{ pid }}.unique()).set{mafft_input}

process mafft_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { 'mafft' }

    publishDir "results/alignment/mafft_{{ pid }}/"

    input:
    file(assembly) from mafft_input.collect()

    output:
    file ("*.align") into {{ output_channel }}
    {% with task_name="mafft", sample_id="val('single')" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    cat ${assembly} > all_assemblies.fasta

    mafft --adjustdirection --thread $task.cpus --auto all_assemblies.fasta > ${workflow.scriptName}.align
    """

}

{{ forks }}