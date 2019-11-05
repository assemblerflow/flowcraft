process pneumocat_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "results/typing/pneumocat_{{ pid }}/$sample_id"

    input:
    set sample_id, file(fastq) from {{ input_channel }}

    output:
    set sample_id, file("*.results.xml") into LOG_pneumocat_{{ pid }}
    file("*.txt")
    {% with task_name="pneumocat" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        #rename input files for pneumocat
        mv ${fastq[0]} ${sample_id}_1.fastq.gz
        mv ${fastq[1]} ${sample_id}_2.fastq.gz

        # run pneumocat
        python /NGStools/PneumoCaT/PneumoCaT.py --fastq_1 ${sample_id}_1.fastq.gz  --fastq_2 ${sample_id}_2.fastq.gz  --threads $task.cpus --cleanup
        mv pneumo_capsular_typing/* . && rm pneumo_capsular_typing/
        echo pass > .status

    } || {
        echo fail > .status
    }
    """

}

process report_pneumocat_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    input:
    set sample_id, file(res) from LOG_pneumocat_{{ pid }}

    output:
    {% with task_name="report_pneumocat" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_pneumocat.py"
}