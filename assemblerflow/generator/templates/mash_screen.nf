
// check if noWinner is provided or not
winnerVar = (params.noWinner == false) ? "-w" : ""

// process to run mashScreen and sort the output into
// sortedMashScreenResults_{sampleId}.txt
process mashScreen_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { "running mash screen for sample: " + sample_id }

    input:
    set sample_id, file(reads) from {{ input_channel }}
    val refFile from IN_reference_file

    output:
    set sample_id, file("sortedMashScreenResults_${id}.txt") into mashScreenResults_{{ pid }}
    {% with task_name="mashScreen", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    mash screen -i ${params.identity} -v ${params.pValue} -p \
    ${task.cpus} ${winnerVar} ${refFile} ${reads} > mashScreenResults_${sample_id}.txt
    sort -gr mashScreenResults_${sample_id}.txt > sortedMashScreenResults_${sample_id}.txt
    """
}

// process to parse the output to json format
process mashOutputJson_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { "dumping json file from: " + mashtxt }

    publishDir 'results/mashscreen/'

    input:
    set sample_id, file(mashtxt) from mashScreenResults_{{ pid }}

    output:
    set sample_id, file("sortedMashScreenResults_${sample_id}.json") optional true into mashScreenOutputChannel_{{ pid }}
    {% with task_name="mashOutputJson", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "mashscreen2json.py"
}

{{ forks }}