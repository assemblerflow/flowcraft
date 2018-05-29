
// runs metaprob
process metaProb_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir 'results/metaprob/'

    input:
    set sample_id, file(reads) from {{ input_channel }}

    output:
    set sample_id, file("*clusters.csv") into metaProbOutChannel_{{ pid }}
    {% with task_name="metaProb", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    MetaProb -pi ${reads} -feature ${params.feature} -m ${params.metaProbQMer}
    """

}