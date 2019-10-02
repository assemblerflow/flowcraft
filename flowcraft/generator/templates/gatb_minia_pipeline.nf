IN_kmers_step_{{ pid }} = Channel.value(params.step{{ param_id }})

IN_besst_iter_{{ pid }} = Channel.value(params.besst_iter{{ param_id }})

error_correction = params.no_error_correction{{ param_id }} ? "true" : "false"
IN_error_correction_{{ pid }} = Channel.value(error_correction)

process GATBMiniaPipeline_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/GATBMiniaPipeline_{{ pid }}/', pattern: '*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val kmer_step from IN_kmers_step_{{ pid }}
    val skip_error_correction from IN_error_correction_{{ pid }}
    val besst_iter from IN_besst_iter_{{ pid }}

    output:
    set sample_id, file('*.fasta') into {{ output_channel }}
    {% with task_name="GATBMiniaPipeline" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o ${sample_id}_GATBMiniaPipeline \
            --step ${kmer_step} --nb-cores $task.cpus --besst_iter ${besst_iter}

        if [ -f "${sample_id}_GATBMiniaPipeline.fasta" ]; then
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