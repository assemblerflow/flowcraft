process metaphlan_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/metaphlan/", pattern: "*.txt"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}

    output:
     set sample_id, file("${sample_id}_profiled_metagenome.txt") into {{ output_channel }}
    {% with task_name="metaphlan" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    metaphlan2.py ${fastq_pair[0]},${fastq_pair[1]} --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq > profiled_metagenome.txt

    mv profiled_metagenome.txt ${sample_id}_profiled_metagenome.txt
    """
}

{{ forks }}