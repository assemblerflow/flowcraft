if ( !params.schemaPath{{ param_id }} ){
    exit 1, "'schemaPath{{ param_id }}' parameter missing"
}
if ( params.chewbbacaTraining{{ param_id }}){
    if (!file(params.chewbbacaTraining{{ param_id }}).exists()) {
        exit 1, "'chewbbacaTraining{{ param_id }}' file was not found: '${params.chewbbacaTraining{{ param_id }}}'"
    }
}
if ( params.schemaSelectedLoci{{ param_id }}){
    if (!file(params.schemaSelectedLoci{{ param_id }}).exists()) {
        exit 1, "'schemaSelectedLoci{{ param_id }}' file was not found: '${params.schemaSelectedLoci{{ param_id }}}'"
    }
}
if ( params.schemaCore{{ param_id }}){
    if (!file(params.schemaCore{{ param_id }}).exists()) {
        exit 1, "'schemaCore{{ param_id }}' file was not found: '${params.schemaCore{{ param_id }}}'"
    }
}

IN_schema_{{ pid }} = Channel.fromPath(params.schemaPath{{ param_id }})


if (params.chewbbacaJson{{ param_id }} == true){
    jsonOpt = "--json"
} else {
    jsonOpt = ""
}

if (params.chewbbacaTraining{{ param_id }}){
    training = "--ptf ${params.chewbbacaTraining{{ param_id }}}"
} else {
    training = ""
}

// If chewbbaca is executed in batch mode, wait for all assembly files
// to be collected on the input channel, and only then execute chewbbaca
// providing all samples simultaneously
if (params.chewbbacaBatch{{ param_id }}) {
    process chewbbaca_batch_{{ pid }} {

        {% include "post.txt" ignore missing %}
        maxForks 1
        scratch false
        if (params.chewbbacaQueue{{ param_id }} != null) {
            queue "${params.chewbbacaQueue{{ param_id}}}"
        }
        publishDir "results/chewbbaca_alleleCall_{{ pid }}/", mode: "copy"

        input:
        file assembly from {{ input_channel }}.map{ it[1] }.collect()
        each file(schema) from IN_schema_{{ pid }}

        output:
        file 'chew_results*'
        file 'cgMLST.tsv' optional true into chewbbacaProfile_{{ pid }}
        {% with task_name="chewbbaca", sample_id="val('single')" %}
        {%- include "compiler_channels.txt" ignore missing -%}
        {% endwith %}

        script:
        """
        {
            set -x
            if [ -d "$schema/temp" ];
            then
                rm -r $schema/temp
            fi

            if [ "$params.schemaSelectedLoci{{ param_id }}" = "null" ];
            then
                inputGenomes=$schema
            else
                inputGenomes=${params.schemaSelectedLoci{{ param_id }}}
            fi

            echo $assembly | tr " " "\n" >> input_file.txt
            chewBBACA.py AlleleCall -i input_file.txt -g \$inputGenomes -o chew_results $jsonOpt --cpu $task.cpus $training
            if [ "$jsonOpt" = "--json" ]; then
                merge_json.py ${params.schemaCore{{ param_id }}} chew_results_*/*/results*
            else
                cp chew_results*/*/results_alleles.tsv cgMLST.tsv
            fi
        } || {
            echo fail > .status
        }
        """
    }

} else {
    process chewbbaca_{{ pid }} {

        // Send POST request to platform
        {% include "post.txt" ignore missing %}

        maxForks 1
        tag { sample_id }
        scratch true
        if (params.chewbbacaQueue{{ param_id }} != null) {
            queue "${params.chewbbacaQueue{{ param_id }}"
        }
        publishDir "results/chewbbaca_alleleCall_{{ pid }}/", mode: "copy"

        input:
        set sample_id, file(assembly) from {{ input_channel }}
        each file(schema) from IN_schema_{{ pid }}

        output:
        file 'chew_results*'
        file '*_cgMLST.tsv' optional true into chewbbacaProfile_{{ pid }}
        {% with task_name="chewbbaca" %}
        {%- include "compiler_channels.txt" ignore missing -%}
        {% endwith %}

        script:
        """
        {
            set -x
            if [ -d "$schema/temp" ];
            then
                rm -r $schema/temp
            fi

            if [ "$params.schemaSelectedLoci{{ param_id }}" = "null" ];
            then
                inputGenomes=$schema
            else
                inputGenomes=${params.schemaSelectedLoci{{ param_id }}}
            fi

            echo $assembly >> input_file.txt
            chewBBACA.py AlleleCall -i input_file.txt -g \$inputGenomes -o chew_results_${sample_id} $jsonOpt --cpu $task.cpus $training --fc
            if [ "$jsonOpt" = "--json" ]; then
                merge_json.py ${params.schemaCore{{ param_id }}} chew_results_*/*/results*
            else
                mv chew_results_*/*/results_alleles.tsv ${sample_id}_cgMLST.tsv
            fi
        } || {
            echo fail > .status
        }
        """
    }
}


process chewbbacaExtractMLST_{{ pid }} {

    publishDir "results/chewbbaca_{{ pid }}/", mode: "copy", overwrite: true

    input:
    file profiles from chewbbacaProfile_{{ pid }}.collect()

    output:
    file "results/cgMLST.tsv"

    """
    head -n1 ${profiles[0]} > chewbbaca_profiles.tsv
    awk 'FNR == 2' $profiles >> chewbbaca_profiles.tsv
    chewBBACA.py ExtractCgMLST -i chewbbaca_profiles.tsv -o results -p $params.chewbbacaProfilePercentage{{ param_id }}
    """

}
