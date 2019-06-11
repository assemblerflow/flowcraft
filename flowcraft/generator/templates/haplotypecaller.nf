haplotypecallerIndexId_{{ pid }} = Channel.value(params.reference{{ param_id }}.split("/").last())
haplotypecallerRef_{{ pid }} = Channel.fromPath("${params.reference{{ param_id }}}.*").collect().toList()
interval_{{ pid }} = Channel.fromPath(params.intervals{{ param_id }})
           .ifEmpty { exit 1, "Interval list file for HaplotypeCaller not found: ${params.intervals}" }
           .splitText()
           .map { it -> it.trim() }

process haplotypecaller_{{ pid }} {

    tag "$interval"

    {% include "post.txt" ignore missing %}

    publishDir "results/variant_calling/haplotypecaller_{{ pid }}"

    input:
    set sample_id, file(bam), file(bai) from {{ input_channel }}
    each interval from interval_{{pid}}
    each file(ref_files) from haplotypecallerRef_{{pid}}
    each index from haplotypecallerIndexId_{{pid}}
   
    output:
    file("*.g.vcf") into haplotypecaller_gvcf
    file("*.g.vcf.idx") into index

    {% with task_name="haplotypecaller", suffix="_${interval}" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    gatk HaplotypeCaller \
      --java-options -Xmx${task.memory.toMega()}M \
      -R ${index}.fasta \
      -O ${sample_id}.g.vcf \
      -I $bam \
      -ERC GVCF \
      -L $interval
    """
}

{{ forks }}