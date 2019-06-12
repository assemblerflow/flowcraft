process mark_duplicates_{{ pid }} {

    {% include "post.txt" ignore missing %}

    input:
    set sample_id, file(bam), file(bai) from {{ input_channel }}
   
    output:
    set val(sample_id), file("${sample_id}_mark_dup.bam"), file("${sample_id}_mark_dup.bai") into {{ output_channel }}
    set file("metrics.txt") into markDupMultiQC_{{pid}}
    {% with task_name="mark_duplicates" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    gatk MarkDuplicates \
      -I $bam \
      -M metrics.txt \
      -O ${sample_id}_mark_dup.bam \
      --CREATE_INDEX
    """
}

{{ forks }}