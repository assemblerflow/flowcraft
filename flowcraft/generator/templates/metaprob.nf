IN_feature_{{ pid }} = Channel.value(params.feature{{ param_id }})
IN_metaProbQMer_{{ pid }} = Channel.value(params.metaProbQMer{{ param_id }})

// runs metaProb
process metaProb_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir 'results/metaprob/'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val feature from IN_feature_{{ pid }}
    val metaProbQMer from IN_metaProbQMer_{{ pid }}

    output:
    set sample_id, file("*clusters.csv") into metaProbOutChannel_{{ pid }}
    {% with task_name="metaProb", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    gunzip -c ${fastq_pair[0]} > ${sample_id}_read1.fastq
    gunzip -c ${fastq_pair[1]} > ${sample_id}_read2.fastq

    MetaProb -pi ${sample_id}_read1.fastq ${sample_id}_read2.fastq -feature ${feature} -m ${pmetaProbQMer}
    """

}

{{ forks }}