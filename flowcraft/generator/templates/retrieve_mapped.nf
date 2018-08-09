process retrieve_mapped_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/mapping/retrieve_mapped_{{ pid }}/'

    input:
    set sample_id, file(bam) from {{ input_channel }}

    output:
    set sample_id , file("*.headersRenamed_*.fq.gz") into {{ output_channel }}
    {% with task_name="retrieve_mapped" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    samtools view -buh -F 12 -o ${sample_id}_samtools.bam -@ $task.cpus ${bam}

    samtools fastq -1 ${sample_id}_mapped_1.fq -2 ${sample_id}_mapped_2.fq ${sample_id}_samtools.bam

    renamePE_samtoolsFASTQ.py -1 ${sample_id}_mapped_1.fq -2 ${sample_id}_mapped_2.fq

    gzip *.headersRenamed_*.fq
    """
}

{{ forks }}