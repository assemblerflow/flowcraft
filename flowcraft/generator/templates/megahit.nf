if ( params.megahitKmers{{ param_id }}.toString().split(" ").size() <= 1 ){
    if (params.megahitKmers{{ param_id }}.toString() != 'auto'){
        exit 1, "'megahitKmers{{ param_id }}' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.megahitKmers{{ param_id }}}"
    }
}
IN_megahit_kmers_{{ pid }} = Channel.value(params.megahitKmers{{ param_id }})

clear = params.clearInput{{ param_id }} ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process megahit_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/megahit_{{ pid }}/', pattern: '*_megahit*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from {{ input_channel }}.join(SIDE_max_len_{{ pid }})
    val kmers from IN_megahit_kmers_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, file('*megahit*.fasta') into {{ output_channel }}
    set sample_id, file('megahit/intermediate_contigs/k*.contigs.fa') into IN_fastg{{ pid }}
    {% with task_name="megahit" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "megahit.py"

}

fastg = params.fastg{{ param_id }} ? "true" : "false"
process megahit_fastg_{{ pid }}{

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "results/assembly/megahit_{{ pid }}/$sample_id", pattern: "*.fastg"

    input:
    set sample_id, file(kmer_files) from IN_fastg{{ pid }}
    val run_fastg from fastg

    output:
    file "*.fastg" optional true
    {% with task_name="megahit_fastg" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    if [ ${run_fastg} == "true" ]
    then
        for kmer_file in ${kmer_files};
        do
            echo \$kmer_file
            k=\$(echo \$kmer_file | cut -d '.' -f 1);
            echo \$k
            megahit_toolkit contig2fastg \$k \$kmer_file > \$kmer_file'.fastg';
        done
    fi
    """
}

{{ forks }}