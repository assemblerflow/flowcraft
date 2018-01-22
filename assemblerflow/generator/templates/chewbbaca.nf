
process chewbbaca {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    maxForks 1
    tag { fastq_id + " getStats" }
    scratch true
    publishDir "results/chewbbaca/${fastq_id}"
    if (params.chewbbacaQueue != null) {
        queue '${params.chewbbacaQueue}'
    }

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    each file(schema) from Channel.fromPath(params.schemaPath)

    output:
    file 'chew_results'
    set fastq_id, val("chewbbaca"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}
    file '.report.json'

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
        chewBBACA.py AlleleCall -i input_file.txt -g ${params.schemaSelectedLoci} -o chew_results --json --cpu $task.cpus -t "${params.chewbbacaSpecies}"
        merge_json.py ${params.schemaCore} chew_results/*/results*
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}

