if (params.index == null && params.reference == null){
    exit 1, "An index or a reference fasta file must be provided."
} else if (params.index != null && params.reference != null){
    exit 1, "Provide only an index OR a reference fasta file."
}

if (params.index){
    index = params.index
} else {
    index = "ref_index"
}

if (params.reference){

    reference_in = Channel.fromPath(params.reference)
        .map{it -> file(it).exists() ? [it.toString().tokenize('/').last().tokenize('.')[0..-2].join('.') ,it] : null}
        .ifEmpty{ exit 1, "No fasta file was provided"}

    process bowtie_build_{{ pid }} {

        // Send POST request to platform
        {% include "post.txt" ignore missing %}

        tag { sample_id }
        storeDir 'index/'
        maxForks 1

        input:
        set sample_id, file(fasta) from reference_in

        output:
        file "ref_index*" into bowtiebuild_index

        script:
        """
        bowtie2-build ${fasta} ref_index > ${sample_id}_bowtie2_build.log
        """
    }
} else {
    bowtiebuild_index = Channel.value(params.index)
}


process bowtie_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
     publishDir 'results/mapping/bowtie_{{ pid }}/'

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    file bowtie2Index from bowtiebuild_index

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