haplotypecallerIndexId_{{ pid }} = Channel.value(params.reference{{ param_id }}.split("/").last())
haplotypecallerRef_{{ pid }} = Channel.fromPath("${params.reference{{ param_id }}}.*").collect().toList()
interval_{{ pid }} = Channel.fromPath(params.intervals{{ param_id }})
           .ifEmpty { exit 1, "Interval list file for HaplotypeCaller not found: ${params.intervals}" }
           .splitText()
           .map { it -> it.trim() }

process haplotypecaller_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag "$interval"

    input:
    set sample_id, file(bam), file(bai) from {{ input_channel }}
    each interval from interval_{{pid}}
    each file(ref_files) from haplotypecallerRef_{{pid}}
    each index from haplotypecallerIndexId_{{pid}}
   
    output:
    file("*.vcf") into haplotypecallerGvcf
    file("*.vcf.idx") into gvcfIndex
    val(sample_id) into sampleId

    {% with task_name="haplotypecaller", suffix="_${interval}" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    gatk HaplotypeCaller \
      --java-options -Xmx${task.memory.toMega()}M \
      -R ${index}.fasta \
      -O ${sample_id}.vcf \
      -I $bam \
      -L $interval
    """
}

process merge_vcfs_{{ pid }} {

    {% include "post.txt" ignore missing %}

    publishDir "results/variant_calling/merge_vcfs_{{ pid }}"

    tag { sample_id }

    input:
    file('*.vcf') from haplotypecallerGvcf.collect()
    file('*.vcf.idx') from gvcfIndex.collect()
    val(sample_id) from sampleId.first()

    output:
    set file("${sample_id}.vcf.gz"), file("${sample_id}.vcf.gz.tbi") into {{ output_channel }}
    {% with task_name="merge_vcfs" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    ## make list of input variant files
    for vcf in \$(ls *vcf); do
      echo \$vcf >> input_variant_files.list
    done

    gatk MergeVcfs \
      --INPUT= input_variant_files.list \
      --OUTPUT= ${sample_id}.vcf.gz
    """

}

{{ forks }}