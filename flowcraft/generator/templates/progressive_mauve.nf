process progressive_mauve_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { 'progressive_mauve' }

    publishDir "results/alignment/progressive_mauve_{{ pid }}/", pattern: '*.align*', mode: 'copy'

    input:
    file(assembly) from {{ input_channel }}.map{ it[1] }.collect()

    output:
    file ("*.align") into {{ output_channel }}
    {% with task_name="progressive_mauve", sample_id="val('single')" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    progressiveMauve --output=${workflow.scriptName}.align --collinear ${assembly}
    """

}

{{ forks }}