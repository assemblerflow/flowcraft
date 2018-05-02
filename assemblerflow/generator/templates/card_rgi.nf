process card_rgi_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/card_rgi/", pattern: "*.txt"

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    val alignmetTool from IN_alignment_tool

    output:
    file("${sample_id}_card_rgi.txt")
    {% with task_name="card_rgi" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    rgi main --input_sequence ${assembly} --output_file ${sample_id}_card_rgi --input_type contig --alignment_tool ${alignmetTool} --low_quality --include_loose -d wgs --clean
    """
}

{{ forks }}