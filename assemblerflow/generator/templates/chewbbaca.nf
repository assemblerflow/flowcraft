
if (params.chewbbacaJson == true){
    jsonOpt = "--json"
} else {
    jsonOpt = ""
}

if (params.chewbbacaTraining){
    training = "--ptf ${params.chewbbacaTraining}"
} else {
    training = ""
}

// If chewbbaca is executed in batch mode, wait for all assembly files
// to be collected on the input channel, and only then execute chewbbaca
// providing all samples simultaneously
if (params.chewbbacaBatch) {
    process chewbbaca_batch_{{ pid }} {

        {% include "post.txt" ignore missing %}
        maxForks 1
        scratch false
        if (params.chewbbacaQueue != null) {
            queue "${params.chewbbacaQueue}"
        }
        publishDir "results/chewbbaca_alleleCall_{{ pid }}/", mode: "copy"

        input:
        file assembly from {{ input_channel }}.map{ it[1] }.collect()
        each file(schema) from IN_schema

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

            if [ "$params.schemaSelectedLoci" = "null" ];
            then
                inputGenomes=$schema
            else
                inputGenomes=${params.schemaSelectedLoci}
            fi

            echo $assembly | tr " " "\n" >> input_file.txt
            chewBBACA.py AlleleCall -i input_file.txt -g \$inputGenomes -o chew_results $jsonOpt --cpu $task.cpus $training
            if [ "$jsonOpt" = "--json" ]; then
                merge_json.py ${params.schemaCore} chew_results_*/*/results*
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
        if (params.chewbbacaQueue != null) {
            queue "${params.chewbbacaQueue}"
        }
        publishDir "results/chewbbaca_alleleCall_{{ pid }}/", mode: "copy"

        input:
        set sample_id, file(assembly) from {{ input_channel }}
        each file(schema) from IN_schema

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

            if [ "$params.schemaSelectedLoci" = "null" ];
            then
                inputGenomes=$schema
            else
                inputGenomes=${params.schemaSelectedLoci}
            fi

            echo $assembly >> input_file.txt
            chewBBACA.py AlleleCall -i input_file.txt -g \$inputGenomes -o chew_results_${sample_id} $jsonOpt --cpu $task.cpus $training --fc
            if [ "$jsonOpt" = "--json" ]; then
                merge_json.py ${params.schemaCore} chew_results_*/*/results*
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
    chewBBACA.py ExtractCgMLST -i chewbbaca_profiles.tsv -o results -p $params.chewbbacaProfilePercentage
    """

}
