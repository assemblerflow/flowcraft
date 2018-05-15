process midas_species_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/annotation/midas/", pattern: "*.txt"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val midasDB from IN_midas_DB

    output:
    file("${sample_id}_midas_species_profile.txt")
    {% with task_name="midas_species" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    run_midas.py species midas/ -d ${midasDB} -t $task.cpus -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --remove_temp

    mv midas/species_profile.txt ./${sample_id}_midas_species_profile.txt
    """
}

{{ forks }}