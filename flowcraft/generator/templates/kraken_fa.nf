IN_kraken_DB_{{ pid }} = Channel.value(params.krakenDB{{ param_id }})


//Process to run Kraken
process kraken_fa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/kraken/", pattern: "*.txt"

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    val krakenDB from IN_kraken_DB_{{ pid }}

    output:
    set sample_id, file("${sample_id}_kraken_report.txt") into toProcess_{{ pid }}
    {% with task_name="kraken_fa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    kraken --preload --db ${krakenDB} --threads $task.cpus --output ${sample_id}_kraken.txt ${assembly}

    kraken-report --db ${krakenDB} ${sample_id}_kraken.txt > ${sample_id}_kraken_report.txt
    """
}


process process_kraken_fa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }


    input:
    set sample_id, file(kraken_report) from toProcess_{{ pid }}

    output:
    //set sample_id, file("${sample_id}_kraken.tsv") into toReport_{{ pid }}
    {% with task_name="process_kraken_fa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_kraken.py"

}

{{ forks }}

