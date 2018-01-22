
process process_spades {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    // This process can only use a single CPU
    cpus 1
    publishDir "reports/assembly/spades_filter", pattern: '*.report.csv', mode: 'copy'

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    val opts from IN_process_spades_opts
    val gsize from IN_genome_size

    output:
    set fastq_id, file('*.assembly.fasta') optional true into {{ output_channel }}
    set fastq_id, val("process_spades"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}
    file '*.report.csv' optional true
    file ".report.json"

    when:
    params.stopAt != "process_spades"

    script:
    template "process_spades.py"

}

{{ forks }}

