
process prokka_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    publishDir "results/annotation/prokka_{{ pid }}/${fastq_id}"

    input:
    set fastq_id, file(assembly) from {{ input_channel }}

    output:
    file "${fastq_id}/*"
    {% with task_name="prokka" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        prokka --outdir $fastq_id --cpus $task.cpus --centre UMMI --compliant \
               --increment 10 $assembly >> .command.log 2>&1
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


