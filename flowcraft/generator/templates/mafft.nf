process mafft_{{ pid }} {

    {% include "post.txt" ignore missing %}

    publishDir "results/mafft/"

    input:
    file(assembly) from {{ input_channel }}.map{ it[1] }.collect()

    output:
    file ("out.align") into {{ output_channel }}
    {% with task_name="mafft", sample_id="val('single')" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """

    cat ${assembly} > all_assemblies.fasta

    mafft --thread $task.cpus --auto all_assemblies.fasta > out.align
    """

}

{{ forks }}