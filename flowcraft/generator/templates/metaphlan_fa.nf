process metaphlan_fa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/metaphlan/", pattern: "*.txt"

    input:
    set sample_id, file(fasta) from {{ input_channel }}

    output:
    set sample_id, file("${sample_id}_krona.txt") into {{ output_channel }}
    file("${sample_id}_profiled_metagenome.txt")
    {% with task_name="metaphlan_fa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    metaphlan2.py ${fasta} --nproc $task.cpus --input_type fasta > ${sample_id}_profiled_metagenome.txt

    metaphlan2krona.py -p ${sample_id}_profiled_metagenome.txt -k ${sample_id}_krona.txt

    parse_krona.py ${sample_id}_krona.txt

    json_str="{'kronaPlot':[{'sample':'${sample_id}','value':\$(cat parsed_krona.txt)}]}"
    echo \$json_str > .report.json
    """
}

{{ forks }}