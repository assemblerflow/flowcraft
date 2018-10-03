process metaphlan_fq_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/metaphlan/", pattern: "*.txt"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}

    output:
     set sample_id, file("${sample_id}_profiled_metagenome.txt") into {{ output_channel }}
    {% with task_name="metaphlan_fq" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    metaphlan2.py ${fastq_pair[0]},${fastq_pair[1]} --bowtie2out metagenome.bowtie2.bz2 --nproc $task.cpus --input_type fastq > ${sample_id}_profiled_metagenome.txt
    """
}

{{ forks }}