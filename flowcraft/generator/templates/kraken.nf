IN_kraken_DB_{{ pid }} = Channel.value(params.krakenDB{{ param_id }})


//Process to run Kraken
process kraken_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/taxonomy/kraken/", pattern: "*.txt"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val krakenDB from IN_kraken_DB_{{ pid }}

    output:
    file("${sample_id}_kraken_report.txt")
    {% with task_name="kraken" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    kraken --preload --fastq-input --db ${krakenDB} --threads $task.cpus --output ${sample_id}_kraken.txt --paired --gzip-compressed ${fastq_pair[0]} ${fastq_pair[1]}

    kraken-report --db ${krakenDB} ${sample_id}_kraken.txt > ${sample_id}_kraken_report.txt
    """
}

{{ forks }}