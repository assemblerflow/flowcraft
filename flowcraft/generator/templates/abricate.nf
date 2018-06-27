
process abricate_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { "${sample_id} ${db}" }
    publishDir "results/annotation/abricate_{{ pid }}/${sample_id}"

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    each db from params.abricateDatabases{{ param_id }}

    output:
    file '*.tsv' into abricate_out_{{ pid }}
    {% with task_name="abricate", suffix="_$db" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # Run abricate
        abricate --db $db $assembly > ${sample_id}_abr_${db}.tsv
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


process process_abricate_{{ pid }} {

    tag "process_abricate_{{ pid }}"

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "report_post.txt" ignore missing %}
    {% endwith %}

    input:
    file abricate_file from abricate_out_{{ pid }}.collect()

    output:
    {% with task_name="process_abricate", sample_id="val('process_abricate')" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_abricate.py"


}



