
process momps_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/typing/momps_{{ pid }}/'

    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})

    output:
    {% with task_name="momps" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    momps.pl -r ${fastq[0]} -f ${fastq[1]} -a $assembly -o res -p $sample_id
    """

}

{{ forks }}

