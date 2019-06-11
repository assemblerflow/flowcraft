process mark_duplicates_{{ pid }} {

    {% include "post.txt" ignore missing %}

    publishDir "results/mapping/mark_duplicates_{{ pid }}"

    input:
    set sample_id, file(bam), file(bai) from {{ input_channel }}
   
    output:
    set sample_id, file("${sample_id}_mark_dup.bam") into {{ output_channel }}
    set file("metrics.txt") into markDupMultiQC_{{pid}}
    {% with task_name="mark_duplicates" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    gatk MarkDuplicates \
      -I $bam \
      -M metrics.txt \
      -O ${sample_id}_mark_dup.bam
    """
}

{{ forks }}