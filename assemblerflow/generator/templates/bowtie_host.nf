process bowtie_host_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from {{ input_channel }}
    val bowtie2Index from IN_index_files

    output:
    set fastq_id , "${fastq_id}*.headersRenamed_*.fq.gz" into {{ output_channel }}
    {% with task_name="bowtie_host" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    bowtie2 -x ${bowtie2Index} -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p task.cpus > ${fastq_id}.bam

    samtools view -buh -f 12 -o ${fastq_id}_samtools.bam -@ 2 ${fastq_id}.bam

    samtools fastq -1 ${fastq_id}_unmapped_1.fq -2 ${fastq_id}_unmapped_2.fq ${fastq_id}_samtools.bam

    renamePE_samtoolsFASTQ.py -1 ${fastq_id}_unmapped_1.fq -2 ${fastq_id}_unmapped_2.fq

    gzip *.headersRenamed_*.fq
    """
}

{{ forks }}