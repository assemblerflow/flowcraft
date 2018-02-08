
process skesa {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }
    publishDir 'results/assembly/skesa', pattern: '*_skesa.assembly.fasta', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair) from {{ input_channel }}

    output:
    set fastq_id, file('*_skesa.assembly.fasta') optional true into {{ output_channel }}
    set fastq_id, val("skesa"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}

    script:
    """
    {
        skesa --fastq ${fastq_pair[0]},${fastq_pair[1]} --gz --use_paired_ends --cores ${task.cpus} > ${fastq_id}_skesa.assembly.fasta
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}