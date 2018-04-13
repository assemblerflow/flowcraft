
// process that runs bowtie2
process mappingBowtie {

    {% include "post.txt" ignore missing %}

    tag { "mapping sample: " + id}

    input:
    set id, file(reads) from {{ input_channel }}
    val bowtie2Index from IN_index_files

    output:
    set id, file("mappingBowtie_${id}.sam") into bowtieResults
    {% with task_name="mappingBowtie", sample_id="id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:

    //if (params.singleEnd == true) {
    //    readsString = "-U ${reads}"
    //}
    //else {
    readsString = "-1 ${reads[0]} -2 ${reads[1]}"
    //}

    """
    bowtie2 -x ${bowtie2Index} ${readsString} -p ${task.cpus} -k \
    ${params.max_k} -5 ${params.trim5} -S mappingBowtie_${id}.sam
    """
}

/**
* samtools faidx is escaped because index file is already provided in docker
* image.
*/
process samtoolsView {

    {% include "post.txt" ignore missing %}

    tag { "samtools commands: " +  id }

    input:
    set id, file(samtoolsFile) from bowtieResults
    val samtoolsIdx from IN_samtools_indexes

    output:
    set id, file("samtoolsDepthOutput_${id}.txt") into samtoolsResults
    {% with task_name="samtoolsView", sample_id="id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    samtools view -b -t ${samtoolsIdx} -@ ${task.cpus} ${samtoolsFile} | \
    samtools sort -@ ${task.cpus} -o samtoolsSorted_${id}.bam
    samtools index samtoolsSorted_${id}.bam
    samtools depth samtoolsSorted_${id}.bam > samtoolsDepthOutput_${id}.txt
    """
}

/**
* These dumping process parses the depth file for each sample and filters it
* depending on the cutoff set by the user.
*/
process jsonDumpingMapping {

    {% include "post.txt" ignore missing %}

    tag { "Dumping json: " +  id }

    publishDir 'results/mapping/'

    input:
    set id, file(depthFile) from samtoolsResults
    val lengthJson from IN_length_json

    output:
    set id, file("samtoolsDepthOutput_${id}.txt_mapping.json") optional true into {{ output_channel }}
    {% with task_name="jsonDumpingMapping", sample_id="id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "mapping2json.py"
}
