process velvet_optimiser_{{pid}} {
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/velvet_optimiser_{{pid}}/', pattern: '*fasta'

    input:
    set sample_id, file(fastq_pair) from {{input_channel}}

    output:
    set sample_id, file('*.fasta') into {{output_channel}}
    {% with task_name="velvet_optimiser" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    VelvetOptimiser.pl -v -s $params.hashs{{ param_id}} -e $params.hashe{{ param_id}} -t $task.cpus \
    -f '-shortPaired -fastq.gz -separate ${fastq_pair[0]} ${fastq_pair[1]}'
    """
}

{{forks}}
