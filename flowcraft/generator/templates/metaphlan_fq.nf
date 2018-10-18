process metaphlan_fq_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/taxonomy/metaphlan/", pattern: "*.txt"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}

    output:
    set sample_id, file("${sample_id}_krona.txt") into {{ output_channel }}
    file("${sample_id}_profiled_metagenome.txt")
    set sample_id, file("${sample_id}_krona.txt") into intoReport{{ pid }}
    {% with task_name="metaphlan_fq" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    metaphlan2.py ${fastq_pair[0]},${fastq_pair[1]} --bowtie2out metagenome.bowtie2.bz2 --nproc $task.cpus --input_type fastq > ${sample_id}_profiled_metagenome.txt

    metaphlan2krona.py -p ${sample_id}_profiled_metagenome.txt -k ${sample_id}_krona.txt

    """
}

process report_metaphlan_fq_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(metaphlan) from intoReport{{ pid }}

    output:
    {% with task_name="report_metaphlan_fq" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_metaphlan.py"

}

{{ forks }}