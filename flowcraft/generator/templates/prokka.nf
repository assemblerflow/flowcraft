
IN_centre_{{ pid }} = Channel.value(params.centre{{ param_id }})

IN_kingdom_{{ pid }} = Channel.value(params.kingdom{{ param_id }})

// check if genus is provided or not
genusVar = (params.genus{{ param_id }} == false) ? "" : "--usegenus --genus ${params.genus{{param_id}}} "

process prokka_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "results/annotation/prokka_{{ pid }}/${sample_id}"

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    val centre from IN_centre_{{ pid }}
    val kingdom from IN_kingdom_{{ pid }}

    output:
    file "${sample_id}/*"
    {% with task_name="prokka" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        prokka --outdir $sample_id --cpus $task.cpus --centre ${centre} \
        --compliant --kingdom ${kingdom} ${genusVar} --increment 10 $assembly
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


