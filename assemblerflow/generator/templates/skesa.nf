
process skesa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }
    publishDir 'results/assembly/skesa_{{ pid }}', pattern: '*_skesa.assembly.fasta', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair) from {{ input_channel }}

    output:
    set fastq_id, file('*_skesa.assembly.fasta') into {{ output_channel }}
    {% with task_name="skesa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        skesa --fastq ${fastq_pair[0]},${fastq_pair[1]} --gz --use_paired_ends \
        --cores ${task.cpus} > ${fastq_id}_skesa.assembly.fasta
        # Test if the output assembly file is empty
        if [ -s ${fastq_id}_skesa.assembly.fasta ]
        then
            echo pass > .status
        else
            echo fail > .status
            exit 1
        fi

    } || {
        echo fail > .status
    }
    """

}