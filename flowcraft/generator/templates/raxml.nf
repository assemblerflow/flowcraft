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
    file ("RAxML_bipartitions.*.nf") into into_json_{{ pid }}
    {% with task_name="raxml", sample_id="val('single')" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    raxmlHPC -s ${alignment} -p 12345 -m ${substitution_model} -T $task.cpus -n $workflow.scriptName -f a -x ${seednumber} -N ${bootstrapnumber}

    # Add information to dotfiles
    version_str="[{'program':'raxmlHPC','version':'8.2.11'}]"
    echo \$version_str > .versions
    """

}

process report_raxml_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { 'raxml' }

    input:
    file(newick) from into_json_{{ pid }}

    output:
    {% with task_name="report_raxml", sample_id="val('single')"  %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_newick.py"

}


{{ forks }}