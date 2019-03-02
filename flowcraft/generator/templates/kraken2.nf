IN_kraken2_DB_{{ pid }} = Channel.value(params.kraken2DB{{ param_id }})


//Process to run Kraken2
process kraken2_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/taxonomy/kraken2/", pattern: "*.txt"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val krakenDB from IN_kraken2_DB_{{ pid }}

    output:
    file("${sample_id}_kraken_report.txt")
    {% with task_name="kraken2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    kraken2 --memory-mapping --threads $task.cpus --report ${sample_id}_kraken_report.txt --db ${krakenDB} --paired \
    --gzip-compressed ${fastq_pair[0]} ${fastq_pair[1]}
    """
}

{{ forks }}