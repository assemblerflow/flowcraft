// Check for the presence of absence of both index and fasta reference
if (params.index{{ param_id }} == null && params.reference{{ param_id }} == null){
    exit 1, "An index or a reference fasta file must be provided."
} else if (params.index{{ param_id }} != null && params.reference{{ param_id }} != null){
    exit 1, "Provide only an index OR a reference fasta file."
}


if (params.reference{{ param_id }}){

    reference_in_{{ pid }} = Channel.fromPath(params.reference{{ param_id }})
        .map{it -> file(it).exists() ? [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.') ,it] : null}
        .ifEmpty{ exit 1, "No fasta file was provided"}

    process bowtie_build_{{ pid }} {

        // Send POST request to platform
        {% include "post.txt" ignore missing %}

        tag { build_id }
        storeDir 'bowtie_index/'
        maxForks 1

        input:
        set build_id, file(fasta) from reference_in_{{ pid }}

        output:
        val build_id into bowtieIndexId_{{ pid }}
        file "${build_id}*.bt2" into bowtieIndex_{{ pid }}

        script:
        """
        bowtie2-build ${fasta} $build_id > ${build_id}_bowtie2_build.log
        """
    }
} else {
    bowtieIndexId_{{ pid }} = Channel.value(params.index{{ param_id }}.split("/").last())
    bowtieIndex_{{ pid }} = Channel.fromPath("${params.index{{ param_id }}}*.bt2").collect().toList()
}


process bowtie_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/mapping/bowtie_{{ pid }}/'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    each index from bowtieIndexId_{{pid}}
    each file(index_files) from bowtieIndex_{{ pid }}

    output:
    set sample_id , file("*.bam") into {{ output_channel }}
    file "*_bowtie2.log"
    {% with task_name="bowtie" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    bowtie2 -x $index -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p $task.cpus 1> ${sample_id}.bam 2> ${sample_id}_bowtie2.log
    """
}

{{ forks }}