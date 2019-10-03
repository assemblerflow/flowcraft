process Fq2Fa_Single_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}

    output:
    set sample_id, file('*.fasta') into {{ output_channel }}
    {% with task_name="Fq2Fa_Single" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        reformat.sh in=§{fastq_pair[0]} in2=§{fastq_pair[0]} out=${sample_id}_reads.fasta

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