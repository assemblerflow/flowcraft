bwaIndexId_{{ pid }} = Channel.value(params.bwaIndex{{ param_id }}.split("/").last())
bwaIndex_{{ pid }} = Channel.fromPath("${params.bwaIndex{{ param_id }}}.*").collect().toList()

process bwa_{{ pid }} {

    {% include "post.txt" ignore missing %}

    publishDir "results/mapping/bwa_{{ pid }}"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    each index from bwaIndexId_{{pid}}
    each file(index_file) from bwaIndex_{{pid}}
   
    output:
    set sample_id, file("${sample_id}.bam"), file("${sample_id}.bam.bai") into {{ output_channel }}
    {% with task_name="bwa" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    bwa mem -M -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:Illumina' -t $task.cpus $index $fastq_pair > ${sample_id}.sam
    samtools sort -o ${sample_id}.bam -O BAM ${sample_id}.sam
    samtools index ${sample_id}.bam
    """
}

{{ forks }}