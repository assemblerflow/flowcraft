IN_alignment_tool_{{ pid }} = Channel.value(params.alignmentTool{{ param_id }})


process card_rgi_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/card_rgi/", pattern: "*.txt"

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    val alignmetTool from IN_alignment_tool_{{ pid }}

    output:
    file("${sample_id}_card_rgi.txt")
    {% with task_name="card_rgi" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    # Place card_rgi source in a read/write location for shifter container
    mkdir card_temp && cp -r /usr/local/lib/python3.5/dist-packages/app/ card_temp
    export PYTHONPATH="\$(pwd)/card_temp:\$PATH"

    rgi main --input_sequence ${assembly} --output_file ${sample_id}_card_rgi --input_type contig --alignment_tool ${alignmetTool} --low_quality --include_loose -d wgs --clean
    """
}

{{ forks }}