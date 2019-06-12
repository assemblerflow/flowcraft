baseRecalibratorFasta_{{ pid }} = Channel.value(params.reference{{ param_id }}.split("/").last())
baseRecalibratorRef_{{ pid }} = Channel.fromPath("${params.reference{{ param_id }}}.*").collect().toList()
baseRecalibratorDbsnp_{{ pid }} = Channel.fromPath("${params.dbsnp{{ param_id }}}")
baseRecalibratorDbsnpIdx_{{ pid }} = Channel.fromPath("${params.dbsnpIdx{{ param_id }}}")
baseRecalibratorGoldenIndel_{{ pid }} = Channel.fromPath("${params.goldenIndel{{ param_id }}}")
baseRecalibratorGoldenIndelIdx_{{ pid }} = Channel.fromPath("${params.goldenIndelIdx{{ param_id }}}")

process base_recalibrator_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set val(sample_id), file(bam), file(bai) from {{ input_channel }}
    each file(reference) from baseRecalibratorRef_{{pid}}
    val(fasta) from baseRecalibratorFasta_{{pid}}
    each file(dbsnp) from baseRecalibratorDbsnp_{{pid}}
    each file(dbsnp_idx) from baseRecalibratorDbsnpIdx_{{pid}}
    each file(golden_indel) from baseRecalibratorGoldenIndel_{{pid}}
    each file(golden_indel_idx) from baseRecalibratorGoldenIndelIdx_{{pid}}
    
    output:
    set sample_id, file("${sample_id}_recal_data.table"), file(bam), file(bai) into baserecalibrator_table
    {% with task_name="base_recalibrator" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    # gunzip dbsnp & golden_indel if gzipped
    [[ "\$(file --mime-type $dbsnp | cut -d' ' -f2)" == "application/x-gzip" ]] && gzip -d --force $dbsnp
    dbsnp=\$(basename $dbsnp .gz)
    [[ "\$(file --mime-type $dbsnp_idx | cut -d' ' -f2)" == "application/x-gzip" ]] && gzip -d --force $dbsnp_idx
    [[ "\$(file --mime-type $golden_indel | cut -d' ' -f2)" == "application/x-gzip" ]] && gzip -d --force $golden_indel
    golden_indel=\$(basename $golden_indel .gz)
    [[ "\$(file --mime-type $golden_indel_idx | cut -d' ' -f2)" == "application/x-gzip" ]] && gzip -d --force $golden_indel_idx

    gatk BaseRecalibrator \
      -I $bam \
      --known-sites \$dbsnp \
      --known-sites \$golden_indel \
      -O ${sample_id}_recal_data.table \
      -R ${fasta}.fasta
    """
}


process apply_bqsr_{{ pid }} {

    {% include "post.txt" ignore missing %}

    publishDir "results/mapping/apply_bqsr_{{ pid }}"

    tag { sample_id }

    input:
    set sample_id, file(baserecalibrator_table), file(bam), file(bai) from baserecalibrator_table
    
    output:
    set sample_id, file("${sample_id}_recalibrated.bam"), file("${sample_id}_recalibrated.bai") into {{ output_channel }}
    {% with task_name="apply_bqsr" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    gatk ApplyBQSR \
      -I $bam \
      -bqsr $baserecalibrator_table \
      -O ${sample_id}_recalibrated.bam \
      --create-output-bam-index
    """
}

{{ forks }}