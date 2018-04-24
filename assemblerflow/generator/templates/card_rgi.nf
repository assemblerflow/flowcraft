process card_rgi_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }

    publishDir "results/annotation/card_rgi/", pattern: "*.txt"

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    val alignmetTool from IN_alignment_tool

    output:
    file("${fastq_id}_card_rgi.txt")
    {% with task_name="card_rgi" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    rgi main --input_sequence ${assembly} --output_file ${fastq_id}_card_rgi --input_type contig --alignment_tool ${alignmetTool} --low_quality --include_loose -d wgs --clean
    """
}

{{ forks }}