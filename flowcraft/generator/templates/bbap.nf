process bbap_{{pid}} {
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/bbap_{{pid}}/', pattern: '*.fasta'

    input:
    set sample_id, file(fastq_pair) from {{input_channel}}

    output:
    set sample_id, file('*.fasta') into {{output_channel}}
    {% with task_name="bbap" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    gunzip -c ${fastq_fair[0]} > reads_1.fq
    gunzip -c ${fastq_fair[1]} > reads_2.fq

    perl /NGStools/BBAP/QC_SB_AC_masterPipeline.pl -p /NGStools/BBAP/ -F 1 -o ${sample_id}_BBAP -O . \
    -b /NGStools/ncbi-blast-2.9.0+/ -a $task.cpus ${fastq_pair[0]} -2 ${fastq_pair[1]}
    """"
}

{{forks}}
