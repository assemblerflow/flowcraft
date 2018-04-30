process kraken_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }

    publishDir "results/annotation/krakken/", pattern: "*.txt"

    input:
    set fastq_id, file(fastq_pair) from {{ input_channel }}
    val krakenDB from IN_kraken_DB

    output:
    file("${fastq_id}_kraken_report.txt")
    {% with task_name="kraken" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    kraken --preload --fastq-input --db ${krakenDB} --threads $task.cpus --output ${fastq_id}_kraken.txt --paired --gzip-compressed ${fastq_pair[0]} ${fastq_pair[1]}

    kraken-report --db ${krakenDB} ${fastq_id}_kraken.txt > ${fastq_id}_kraken_report.txt
    """
}

{{ forks }}