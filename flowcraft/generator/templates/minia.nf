IN_kmer_{{ pid }} = Channel.value(params.miniaKmer{{ param_id }})

process Minia_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/Minia_{{ pid }}/', pattern: '*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val kmer from IN_kmers_{{ pid }}

    output:
    set sample_id, file('*.fasta') into {{ output_channel }}
    {% with task_name="Minia" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        cat ${fastq_pair[0]} ${fastq_pair[1]} > ${sample_id}_reads.fq.gz
        minia -out ${sample_id}_minia.fasta -in ${sample_id}_reads.fq.gz -kmer-size ${kmer}

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