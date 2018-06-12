
process prokka_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "results/annotation/prokka_{{ pid }}/${sample_id}"

    input:
    set sample_id, file(assembly) from {{ input_channel }}

    output:
    file "${sample_id}/*"
    {% with task_name="prokka" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        prokka --outdir $sample_id --cpus $task.cpus --centre UMMI --compliant \
               --increment 10 $assembly >> .command.log 2>&1
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


