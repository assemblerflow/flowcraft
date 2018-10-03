process metaphlan_to_krona_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(metaphlan_profile) from {{ input_channel }}

    output:
     set sample_id, file("${sample_id}_krona.txt") into {{ output_channel }}
    {% with task_name="metaphlan_to_krona" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    metaphlan2krona.py -p ${metaphlan_profile} -k ${sample_id}_krona.txt
    """
}

{{ forks }}