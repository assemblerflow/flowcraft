IN_shared_hashes_{{ pid }} = Channel.value(params.shared_hashes{{ param_id }})

if (binding.hasVariable("SIDE_mashSketchOutChannel_{{ pid }}")){
    IN_reference_file_{{ pid }} = SIDE_mashSketchOutChannel_{{ pid }}
} else {
    IN_reference_file_{{ pid }} = Channel.value(params.refFile{{ param_id }})
}

// runs mash dist
process runMashDist_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
//    file refFile from SIDE_mashSketchOutChannel_{{ pid }}.ifEmpty(IN_reference_file_{{ pid }})
    set sample_id, file(fasta) from {{ input_channel }}
    val refFile from IN_reference_file_{{ pid }}

    output:
    set sample_id, fasta, file("*_mashdist.txt") into mashDistOutChannel_{{ pid }}
    {% with task_name="runMashDist", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    mash dist -i -p ${task.cpus} -v ${params.pValue{{ param_id }}} \
    -d ${params.mash_distance{{ param_id }}} ${refFile} ${fasta} > ${fasta}_mashdist.txt
    """

}

// parses mash dist output to a json file that can be imported into pATLAS
process mashDistOutputJson_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir 'results/mashdist/mashdist_json_{{ pid }}/', mode: 'copy'

    input:
    set sample_id, fasta, file(mashtxt) from mashDistOutChannel_{{ pid }}
    val shared_hashes from IN_shared_hashes_{{ pid }}

    output:
    set sample_id, file("*.json") optional true into {{ output_channel }}
    {% with task_name="mashDistOutputJson", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "mashdist2json.py"

}

{{ forks }}
