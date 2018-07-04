
process mlst_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(assembly) from {{ input_channel }}

    output:
    file '*.mlst.txt' into LOG_mlst_{{ pid }}
    set sample_id, file(assembly), file(".status") into MAIN_mlst_out_{{ pid }}
    {% with task_name="mlst" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        expectedSpecies=${params.mlstSpecies{{ param_id }}}
        mlst $assembly >> ${sample_id}.mlst.txt
        mlstSpecies=\$(cat *.mlst.txt | cut -f2)
        json_str="{'expectedSpecies':\'\$expectedSpecies\','species':'\$mlstSpecies','st':'\$(cat *.mlst.txt | cut -f3)','tableRow':[{'sample':'${sample_id}','data':[{'header':'mlst','value':'\$mlstSpecies','table':'typing'}]}]}"
        echo \$json_str > .report.json

        if [ ! \$mlstSpecies = \$expectedSpecies ];
        then
            printf fail > .status
        else
            printf pass > .status
        fi

    } || {
        printf fail > .status
    }
    """
}

process compile_mlst_{{ pid }} {

    publishDir "results/annotation/mlst_{{ pid }}/"

    input:
    file res from LOG_mlst_{{ pid }}.collect()

    output:
    file "mlst_report.tsv"

    script:
    """
    cat $res >> mlst_report.tsv
    """
}

{{ output_channel }} = Channel.create()
MAIN_mlst_out_{{ pid }}
    .filter{ it[2].text != "fail" }
    .map{ [it[0], it[1]] }
    .set{ {{output_channel}} }


{{ forks }}

