process krona_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/krona/", pattern: "*.html"

    input:
    set sample_id, file(krona_file) from {{ input_channel }}

    output:
    set sample_id, file("${sample_id}_krona.html")
    {% with task_name="krona" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    ktImportText ${krona_file} -o ${sample_id}_krona.html

    json_str="{'kronaPlot':[{'sample':'${sample_id}','value':'\$(cat ${sample_id}_krona.html)'}]}"
    echo \$json_str > .report.json

    """
}

{{ forks }}