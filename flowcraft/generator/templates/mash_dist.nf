IN_shared_hashes_{{ pid }} = Channel.value(params.shared_hashes{{ param_id }})

IN_mash_dist_input = Channel.create()
// If the side channel with the sketch exists, join the corresponding .msh file
// with the appropriate sample_id
if (binding.hasVariable("SIDE_mashSketchOutChannel_{{ pid }}")){
    {{ input_channel }}
        .join(SIDE_mashSketchOutChannel_{{ pid }})
        .into(IN_mash_dist_input)
// Otherwise, always use the .msh file provided in the docker image
} else {
    {{ input_channel }}
        .map{ it -> [it[0], it[1], params.refFile{{ param_id }}] }
        .into(IN_mash_dist_input)
}

// runs mash dist
process runMashDist_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir 'results/mashdist/mashdist_txt_{{ pid }}/'

    input:
    set sample_id, file(fasta), refFile from IN_mash_dist_input

    output:
    set sample_id, file(fasta), file("*_mashdist.txt") into mashDistOutChannel_{{ pid }}
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

    publishDir 'results/mashdist/mashdist_json_{{ pid }}/'

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
