IN_substitution_model_{{ pid }} = Channel.value(params.substitutionModel{{ param_id }})


process raxml_{{ pid }} {

    {% include "post.txt" ignore missing %}

    publishDir "results/raxml/"

    input:
    file(alignment) from {{ input_channel }}
    val substitution_model from IN_substitution_model_{{ pid }}

    output:
    file ("*.tree") into {{ output_channel }}
    {% with task_name="raxml", sample_id="val('single')" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    raxmlHPC -s ${alignment} -p 12345 -m ${substitution_model} -T $task.cpus -n tree
    """

}

{{ forks }}