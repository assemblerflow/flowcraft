process maxbin2_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/maxbin2_{{ pid }}/${sample_id}/"

    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})

    output:
    file '*_maxbin.*.fasta' into binCh
    file '*_maxbin.{abundance,log,summary}'
    {% with task_name="maxbin2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    run_MaxBin.pl -contig ${assembly} -out ${sample_id}_maxbin -reads ${fastq[0]} -reads2 ${fastq[1]} -thread $task.cpus -min_contig_length ${params.min_contig_lenght} -max_iteration ${params.max_iteration} -prob_threshold ${params.prob_threshold}
    """
}

{{ output_channel }} = Channel.create()

binCh.flatMap().map{ it -> [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it] }
    .into({{ output_channel }})

{{ forks }}