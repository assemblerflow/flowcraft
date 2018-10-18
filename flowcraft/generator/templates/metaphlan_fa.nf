process metaphlan_fa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/taxonomy/metaphlan/", pattern: "*.txt"

    input:
    set sample_id, file(fasta) from {{ input_channel }}

    output:
    set sample_id, file("${sample_id}_krona.txt") into {{ output_channel }}
    file("${sample_id}_profiled_metagenome.txt")
    set sample_id, file("${sample_id}_krona.txt") into intoReport{{ pid }}
    {% with task_name="metaphlan_fa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    metaphlan2.py ${fasta} --nproc $task.cpus --input_type fasta > ${sample_id}_profiled_metagenome.txt

    metaphlan2krona.py -p ${sample_id}_profiled_metagenome.txt -k ${sample_id}_krona.txt

    """
}

process report_metaphlan_fa_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(metaphlan) from intoReport{{ pid }}

    output:
    {% with task_name="report_metaphlan_fa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_metaphlan.py"

}

{{ forks }}