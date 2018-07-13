if (binding.hasVariable("SIDE_mashSketchOutChannel_{{ pid }}")){
    IN_reference_file_{{ pid }} = SIDE_mashSketchOutChannel_{{ pid }}
} else {
    IN_reference_file_{{ pid }} = Channel.value(params.refFile{{ param_id }})
}

// check if noWinner is provided or not
winnerVar = (params.noWinner{{ param_id }} == false) ? "-w" : ""

// process to run mashScreen and sort the output into
// sortedMashScreenResults_{sampleId}.txt
process mashScreen_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(reads) from {{ input_channel }}
    val refFile from IN_reference_file_{{ pid }}

    output:
    set sample_id, file("sortedMashScreenResults*.txt") into mashScreenResults_{{ pid }}
    {% with task_name="mashScreen", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    mash screen -i ${params.identity{{ param_id }}} -v ${params.pValue{{ param_id }}} -p \
    ${task.cpus} ${winnerVar} ${refFile} ${reads} > mashScreenResults_${sample_id}.txt
    sort -gr mashScreenResults_${sample_id}.txt > sortedMashScreenResults_${sample_id}.txt
    """
}

// process to parse the output to json format
process mashOutputJson_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir 'results/mashscreen/mashscreen_json_{{ pid }}', mode: 'copy'

    input:
    set sample_id, file(mashtxt) from mashScreenResults_{{ pid }}

    output:
    set sample_id, file("sortedMashScreenResults*.json") optional true into mashScreenOutputChannel_{{ pid }}
    {% with task_name="mashOutputJson", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "mashscreen2json.py"
}

{{ forks }}