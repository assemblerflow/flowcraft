IN_min_contig_lenght_{{ pid }} = Channel.value(params.min_contig_lenght{{ param_id }})
IN_max_iteration_{{ pid }} = Channel.value(params.max_iteration{{ param_id }})
IN_prob_threshold_{{ pid }} = Channel.value(params.prob_threshold{{ param_id }})

process maxbin2_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/maxbin2_{{ pid }}/${sample_id}/"

    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})
    val minContigLenght from IN_min_contig_lenght_{{ pid }}
    val maxIterations from IN_max_iteration_{{ pid }}
    val probThreshold from IN_prob_threshold_{{ pid }}


    output:
    file '*_maxbin.*.fasta' into binCh_{{ pid }}
    file '*_maxbin.{abundance,log,summary}'
    {% with task_name="maxbin2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    run_MaxBin.pl -contig ${assembly} -out ${sample_id}_maxbin -reads ${fastq[0]} -reads2 ${fastq[1]} -thread $task.cpus -min_contig_length ${minContigLenght} -max_iteration ${maxIterations} -prob_threshold ${probThreshold}
    """
}

{{ output_channel }} = Channel.create()

binCh_{{ pid }}.flatMap().map{ it -> [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it] }
    .into({{ output_channel }})

{{ forks }}