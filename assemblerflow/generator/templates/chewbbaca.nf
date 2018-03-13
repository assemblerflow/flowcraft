
if (params.chewbbacaToPhyloviz == true){
    jsonOpt = ""
} else {
    jsonOpt = "--json"
}

process chewbbaca {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    maxForks 1
    tag { fastq_id + " getStats" }
    scratch true
    if (params.chewbbacaQueue != null) {
        queue '${params.chewbbacaQueue}'
    }

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    each file(schema) from Channel.fromPath(params.schemaPath)

    output:
    file 'chew_results'
    file '*_cgMLST.tsv' optional true into chewbbacaProfile
    {% with task_name="chewbbaca" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    when:
    params.chewbbacaRun == true

    script:
    """
    {
        if [ -d ${params.schemaPath}/temp ];
        then
            rm -r ${params.schemaPath}/temp
        fi

        echo $assembly >> input_file.txt
        chewBBACA.py AlleleCall -i input_file.txt -g ${params.schemaSelectedLoci} -o chew_results $jsonOpt --cpu $task.cpus -t "${params.chewbbacaSpecies}"
        if [ ! $jsonOpt = ""]; then
            merge_json.py ${params.schemaCore} chew_results/*/results*
        else
            mv chew_results/*/results_alleles.tsv ${fastq_id}_cgMLST.tsv
        fi
    } || {
        echo fail > .status
    }
    """

}


process chewbbacaExtractMLST {

    publishDir "results/chewbbaca/", mode: "copy", overwrite: true

    input:
    file profiles from chewbbacaProfile.collect()

    output:
    file "results/cgMLST.tsv"

    """
    head -n1 ${profiles[0]} > chewbbaca_profiles.tsv
    awk 'FNR == 2' $profiles >> chewbbaca_profiles.tsv
    chewBBACA.py ExtractCgMLST -i chewbbaca_profiles.tsv -o results -p $params.chewbbacaProfilePercentage
    """

}
