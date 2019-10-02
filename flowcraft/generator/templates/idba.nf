
process fastq2fasta_{{pid}}{

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(fastq_pair) from {{input_channel}}

    output:
    set sample_id, file('reads.fasta') into OUT_fastq2fastq_{{pid}}

    script:
    """
    gunzip -c ${fastq_fair[0]} > reads_1.fq
    gunzip -c ${fastq_fair[1]} > reads_2.fq

    fq2fa --merge reads_1.fq reads_2.fq reads.fasta
    """
}

process idba_{{pid}} {
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/idba_{{pid}}/', pattern: '*fasta'

    input:
    set sample_id, file(reads_fasta) from OUT_fastq2fastq_{{pid}}

    output:
    set sample_id, file('*.fasta') into {{output_channel}}
    {% with task_name="idba" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    idba_ud -l ${reads_fasta} --num_threads $task.cpus -o .

    """
}

{{forks}}
