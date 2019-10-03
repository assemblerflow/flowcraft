IN_kmers_{{ pid }} = Channel.value(params.kmers{{ param_id }})

IN_threshold_{{ pid }} = Channel.value(params.threshold{{ param_id }})

IN_algorithm_{{ pid }} = Channel.value(params.algorithm{{ param_id }})

process pandaseq_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/pandaseq_{{ pid }}/', pattern: '*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val kmer from IN_kmers_{{ pid }}
    val threshold_val from IN_threshold_{{ pid }}
    val set_algorithm from IN_algorithm_{{ pid }}

    output:
    set sample_id, file('*.fasta') into {{ output_channel }}
    {% with task_name="pandaseq" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        pandaseq -f ${fastq_pair[0]} -r ${fastq_pair[1]} -B -A ${set_algorithm} -k ${kmer}
        -t ${threshold_val} -T $task.cpus -w ${sample_id}_pandaseq.fasta

        if [ -f "*.fasta" ]; then
            printf pass > .status
        else
            printf fail > .status
        fi

    } || {
     printf fail > .status
    }
    """

}

{{ forks }}