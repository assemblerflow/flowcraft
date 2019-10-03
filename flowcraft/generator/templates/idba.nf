process idba_{{pid}} {
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/idba_{{pid}}/', pattern: '*fasta'

    input:
    set sample_id, file(fasta_reads_single) from {{input_channel}}

    output:
    set sample_id, file('*.fasta') into {{output_channel}}
    {% with task_name="idba" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    idba_ud -l ${fasta_reads_single} --num_threads $task.cpus -o .

    """
}

{{forks}}
