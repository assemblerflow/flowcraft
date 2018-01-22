
process spades {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    publishDir 'results/assembly/spades/', pattern: '*_spades.assembly.fasta', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), max_len from {{ input_channel }}.join(SIDE_max_len_{{ pid }})
    val opts from IN_spades_opts
    val kmers from IN_spades_kmers

    output:
    set fastq_id, file('*_spades.assembly.fasta') optional true into {{ output_channel }}
    set fastq_id, val("spades"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}
    file ".report.json"

    when:
    params.stopAt != "spades"

    script:
    template "spades.py"

}

{{ forks }}

