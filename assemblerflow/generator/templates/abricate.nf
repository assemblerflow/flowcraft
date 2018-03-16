
process abricate_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { "${fastq_id} ${db}" + " getStats"}
    publishDir "results/annotation/abricate_{{ pid }}/${fastq_id}"

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    each db from params.abricateDatabases

    output:
    file '*.tsv' into abricate_out_{{ pid }}
    {% with task_name="abricate" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    when:
    params.abricateRun == true && params.annotationRun

    script:
    """
    {
        # Run abricate
        abricate --db $db $assembly > ${fastq_id}_abr_${db}.tsv
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


process process_abricate_{{ pid }} {

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "report_post.txt" ignore missing %}
    {% endwith %}

    input:
    file abricate_file from abricate_out_{{ pid }}.collect()

    script:
    template "process_abricate.py"


}



