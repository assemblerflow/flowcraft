
IN_genome_size_{{ pid }} = Channel.value(params.genomeSize{{ param_id }})
    .map{it -> it.toString().isNumber() ? it : exit(1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize{{ param_id }}}'")}

IN_depth_{{ pid }} = Channel.value(params.depth{{ param_id }})
    .map{it -> it.toString().isNumber() ? it : exit(1, "The depth parameter must be a number or a float. Provided value: '${params.depth{{ param_id }}}'")}

clear = params.clearAtCheckpoint ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process downsample_fastq_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { "${sample_id}" }
    publishDir "results/downsample_fastq_{{ pid }}/", pattern: "_ss.*"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val gsize from IN_genome_size_{{ pid }}
    val depth from IN_depth_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, file('*_ss.*') into {{ output_channel }}
    {% with task_name="downsample_fastq" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "downsample_fastq.py"

}

{{ forks }}