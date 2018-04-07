
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

process chewbbaca_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    maxForks 1
    tag { fastq_id + " getStats" }
    scratch true
    if (params.chewbbacaQueue != null) {
        queue "${params.chewbbacaQueue}"
    }
    publishDir "results/chewbbaca_alleleCall_{{ pid }}/", mode: "copy"

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    each file(schema) from Channel.fromPath(params.schemaPath)

    output:
    file 'chew_results_*'
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
        chewBBACA.py AlleleCall -i input_file.txt -g \$inputGenomes -o chew_results_${fastq_id} $jsonOpt --cpu $task.cpus $training
        if [ "$jsonOpt" = "--json" ]; then
            merge_json.py ${params.schemaCore} chew_results_*/*/results*
        else
            mv chew_results_*/*/results_alleles.tsv ${fastq_id}_cgMLST.tsv
        fi
    } || {
        echo fail > .status
    }
    """

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
