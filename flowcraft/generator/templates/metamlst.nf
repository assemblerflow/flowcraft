IN_metamlstDB_{{ pid }} = Channel.value(params.metamlstDB{{ param_id }})
IN_metamlstDB_index_{{ pid }} = Channel.value(params.metamlstDB_index{{ param_id }})


process metamlst_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "results/annotation/metamlst_{{ pid }}/${sample_id}", saveAs: { it.split("/").last() }

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val metamlstDB from IN_metamlstDB_{{ pid }}
    val metamlstDB_index from IN_metamlstDB_index_{{ pid }}

    output:
    file 'out/merged/*.txt' optional true
    {% with task_name="metamlst" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    bowtie2 --very-sensitive-local -a --no-unal -x ${metamlstDB_index} -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} | samtools view -bS - > ${sample_id}.bam

    metamlst.py -d ${metamlstDB} ${sample_id}.bam

    metamlst-merge.py -d ${metamlstDB} out/
    """

}

{{ forks }}