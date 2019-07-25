if (params.reference{{ param_id }}) {
    Channel
        .fromPath("${params.reference{{ param_id }}}.fasta")
        .ifEmpty { exit 1, "FASTA annotation file not found: ${params.reference{{ param_id }}}" }
        .set { hisat2Fasta_{{pid}} }
} else if (params.hisat2_index{{ param_id }}) {
    Channel
        .fromPath(params.hisat2_index{{ param_id }})
        .ifEmpty { exit 1, "Folder containing Hisat2 indexes for reference genome not found: ${params.hisat2_index{{ param_id }}}" }
        .set { hisat2Index_{{pid}} }
    hisat2IndexName_{{pid}} = Channel.value( "${params.hisat2_index_name{{ param_id }}}" )
} else {
    exit 1, "Please specify either `--reference /path/to/file_basename` OR `--hisat2_index /path/to/hisat2_index_folder` AND `--hisat2_index_name hisat2_index_folder/basename`"
}

if (!params.hisat2_index{{ param_id }}) {
  process make_hisat2_index_{{ pid }} {

    {% include "post.txt" ignore missing %}
    tag "$fasta"

    input:
    each file(fasta) from hisat2Fasta_{{pid}}

    output:
    val "hisat2_index/${fasta.baseName}.hisat2_index" into hisat2IndexName_{{pid}}
    file "hisat2_index" into hisat2Index_{{pid}}

    """
    mkdir hisat2_index
    hisat2-build -p ${task.cpus} $fasta hisat2_index/${fasta.baseName}.hisat2_index
    """
  }
}

process hisat2_{{ pid }} {

    {% include "post.txt" ignore missing %}
    tag "$sample_id"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    each index_name from hisat2IndexName_{{pid}}
    each file(index) from hisat2Index_{{pid}}

    output:
    set sample_id, file("${sample_id}.sam") into  hisat2Sam_{{pid}}
    {% with task_name="hisat2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    hisat2 \
    -p ${task.cpus} \
    -x $index_name \
    -1 ${fastq_pair[0]} \
    -2 ${fastq_pair[1]} \
    -S ${sample_id}.sam
    """
}

process samtools_sort_{{ pid }} {

    {% include "post.txt" ignore missing %}
    publishDir "results/mapping/hisat2_{{ pid }}"
    tag "$sample_id"

    input:
    set sample_id, file(sam) from hisat2Sam_{{pid}}

    output: 
    set sample_id, file("${sample_id}.sorted.bam"), file("${sample_id}.sorted.bam.bai") into {{ output_channel }}
    {% with task_name="samtools_sort" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    samtools view -Sb $sam > ${sample_id}.bam
    samtools sort -T ${sample_id}.bam.tmp ${sample_id}.bam -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

{{ forks }}