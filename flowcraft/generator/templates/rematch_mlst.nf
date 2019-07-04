if (params.mlstSpecies{{ param_id }} == null){
    exit 1, "A species must be provided."
}

doubleRun = params.doubleRun{{ param_id }} ? "true" : "false"
IN_double_run{{ pid }} = Channel.value(params.doubleRun{{ param_id }})

IN_species_{{ pid }} = Channel.value(params.mlstSpecies{{ param_id }})

process rematch_mlst_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val species from IN_species_{{ pid }}
    val double from IN_double_run{{ pid }}

    output:
    file 'mlst_report_*.tab' into LOG_mlst_{{ pid }}
    set sample_id, file(fastq_pair) into {{output_channel}}
    {% with task_name="rematch_mlst" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # rematch requires the input to be passed like this...
        mkdir reads/
        mkdir reads/${sample_id}
        cp ${fastq_pair[0]} reads/${sample_id}/${sample_id}_1.fq.gz
        cp ${fastq_pair[1]} reads/${sample_id}/${sample_id}_2.fq.gz

        # run rematch
        if [ double = "true" ]
        then
            rematch.py --mlst "${species}" --mlstReference -w reads/ -j $task.cpus --doubleRun --mlstRun second
        else
            rematch.py --mlst "${species}" --mlstReference -w reads/ -j $task.cpus
        fi

        json_str="{'species':'${species}'',\
            'st':'\$(tail -1 */mlst_report* | cut -f4)',\
            'tableRow':[{'sample':'${sample_id}','data':[\
                {'header':'MLST species','value':'${species}','table':'typing'},\
                {'header':'MLST ST','value':'\$(tail -1 */mlst_report* | cut -f4)','table':'typing'}]}]}"
        echo \$json_str > .report.json

        mv */mlst_report* mlst_report_${sample_id}.tab

    } || {
        printf fail > .status
    }
    """
}

process compile_rematch_mlst_{{ pid }} {

    publishDir "results/annotation/rematch_mlst_{{ pid }}/"

    input:
    file res from LOG_mlst_{{ pid }}.collect()

    output:
    file "rematch_mlst_report.tsv"

    script:
    """
    head -1 ${res[0]} > rematch_mlst_report.tsv
    for filename in $res; do
        tail -n 1 "\$filename" >> rematch_mlst_report.tsv
    done

    """
}


{{ forks }}

