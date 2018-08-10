IN_substitution_model_{{ pid }} = Channel.value(params.substitutionModel{{ param_id }})
IN_seed_number_{{ pid }} = Channel.value(params.seedNumber{{ param_id }})
IN_bootstrap_number_{{ pid }} = Channel.value(params.bootstrap{{ param_id }})

process raxml_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { 'raxml' }

    publishDir "results/phylogeny/raxml_{{ pid }}/"

    input:
    file(alignment) from {{ input_channel }}
    val substitution_model from IN_substitution_model_{{ pid }}
    val seednumber from IN_seed_number_{{ pid }}
    val bootstrapnumber from IN_bootstrap_number_{{ pid }}

    output:
    file ("RAxML_*") into {{ output_channel }}
    {% with task_name="raxml", sample_id="val('single')" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    raxmlHPC -s ${alignment} -p 12345 -m ${substitution_model} -T $task.cpus -n $workflow.scriptName -f a -x ${seednumber} -N ${bootstrapnumber}
    """

}

{{ forks }}