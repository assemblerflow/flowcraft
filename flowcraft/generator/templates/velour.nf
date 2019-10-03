IN_kmer_{{ pid }} = Channel.value(params.velourKmer{{ param_id }})

process velour_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/Velour_{{ pid }}/', pattern: 'out/*.fasta', mode: 'copy'

    input:
    set sample_id, file(fasta_reads_pair) from {{ input_channel }}
    val kmer from IN_kmers_{{ pid }}

    output:
    set sample_id, file('out/*.fasta') into {{ output_channel }}
    {% with task_name="velour" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        velour out ${kmer} -file_format fasta -read_type shortPaired ${fasta_reads_pair}

        if [ -f "out/*.fasta" ]; then
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