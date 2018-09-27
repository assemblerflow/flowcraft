IN_index_files_{{ pid }} = Channel.value(params.refIndex{{ param_id }})

clear = params.clearInput{{ param_id }} ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process remove_host_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/mapping/remove_host_{{ pid }}/', pattern: '*_bowtie2.log', mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val bowtie2Index from IN_index_files_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id , file("${sample_id}*.headersRenamed_*.fq.gz") into {{ output_channel }}
    set sample_id, file("*_bowtie2.log") into into_json_{{ pid }}
    {% with task_name="remove_host" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        bowtie2 -x ${bowtie2Index} -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p $task.cpus 1> ${sample_id}.bam 2> ${sample_id}_bowtie2.log

        samtools view -buh -f 12 -o ${sample_id}_samtools.bam -@ $task.cpus ${sample_id}.bam

        rm ${sample_id}.bam

        samtools fastq -1 ${sample_id}_unmapped_1.fq -2 ${sample_id}_unmapped_2.fq ${sample_id}_samtools.bam

        rm ${sample_id}_samtools.bam

        renamePE_samtoolsFASTQ.py -1 ${sample_id}_unmapped_1.fq -2 ${sample_id}_unmapped_2.fq

        gzip *.headersRenamed_*.fq
        rm *.fq

        if [ "$clear" = "true" ];
        then
            work_regex=".*/work/.{2}/.{30}/.*"
            file_source1=\$(readlink -f \$(pwd)/${fastq_pair[0]})
            file_source2=\$(readlink -f \$(pwd)/${fastq_pair[1]})
            if [[ "\$file_source1" =~ \$work_regex ]]; then
                rm \$file_source1 \$file_source2
            fi
        fi

    } || {
        echo fail > .status
    }
    """
}



process report_remove_host_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(bowtie_log) from into_json_{{ pid }}

    output:
    {% with task_name="report_remove_host" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_mapping.py"

}

{{ forks }}