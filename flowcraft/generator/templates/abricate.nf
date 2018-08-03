if ( params.abricateDataDir{{ param_id }} ){
    if ( !file(params.abricateDataDir{{ param_id }}).exists() ){
        exit 1, "'abricateDataDir{{ param_id }}' data directory was not found: '${params.abricateDatabases{{ param_id }}}'"
    }
    dataDirOpt = "--datadir ${params.abricateDataDir{{ param_id }}}"
} else {
    dataDirOpt = ""
}

if ( !params.abricateMinId{{ param_id }}.toString().isNumber() ){
    exit 1, "'abricateMinId{{ param_id }}' parameter must be a number. Provide value: '${params.abricateMinId{{ param_id }}}'"
}

if ( !params.abricateMinCov{{ param_id }}.toString().isNumber() ){
    exit 1, "'abricateMinCov{{ param_id }}' parameter must be a number. Provide value: '${params.abricateMinCov{{ param_id }}}'"
}


process abricate_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { "${sample_id} ${db}" }
    publishDir "results/annotation/abricate_{{ pid }}/${sample_id}"

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    each db from params.abricateDatabases{{ param_id }}
    val min_id from Channel.value(params.abricateMinId{{ param_id }})
    val min_cov from Channel.value(params.abricateMinCov{{ param_id }})

    output:
    file '*.tsv' into abricate_out_{{ pid }}
    {% with task_name="abricate", suffix="_$db" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # Run abricate
        abricate $dataDirOpt --minid $min_id --mincov $min_cov --db $db $assembly > ${sample_id}_abr_${db}.tsv
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



