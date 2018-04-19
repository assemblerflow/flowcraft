
// runs mash dist
process runMashDist_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { "running mash dist for fasta file: " + fasta }

    input:
    set sample_id, file(fasta) from {{ input_channel }}
    val refFile from IN_reference_file

    output:
    set sample_id, fasta, file("${fasta}_mashdist.txt") into mashDistOutChannel_{{ pid }}
    {% with task_name="runMashDist", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    mash dist -i -p ${task.cpus} -v ${params.pValue} \
    -d ${params.mash_distance} ${refFile} ${fasta} > ${fasta}_mashdist.txt
    """

}

// parses mash dist output to a json file that can be imported into pATLAS
process mashDistOutputJson_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { "dumping json file from: " + mashtxt }

    publishDir 'results/mashdist/'

    input:
    set sample_id, fasta, file(mashtxt) from mashDistOutChannel_{{ pid }}

    output:
    set sample_id, file("*.json") optional true into {{ output_channel }}
    {% with task_name="mashDistOutputJson", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "mashdist2json.py"

}

{{ forks }}
