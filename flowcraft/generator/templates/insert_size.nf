IN_plot{{ pid }} = params.distribution_plot{{ param_id }} ? "True" : "False"


process assembly_mapping_statistics_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})

    output:
    set sample_id, file("samtools_stats.txt") into IN_insert_size_{{ pid }}
    {% with task_name="assembly_mapping_statistics" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        echo [DEBUG] BUILDING BOWTIE INDEX FOR ASSEMBLY: $assembly >> .command.log 2>&1
        bowtie2-build --threads ${task.cpus} $assembly genome_index >> .command.log 2>&1

        echo [DEBUG] MAPPING READS FROM $fastq >> .command.log 2>&1
        bowtie2 -q --very-fast --threads ${task.cpus} -x genome_index -1 ${fastq[0]} -2 ${fastq[1]} \
        --fr -I 0 -X 2000 --no-discordant --no-mixed --no-unal -S alignment.sam >> .command.log 2>&1

        echo [DEBUG] GET STATISTICS FROM SAM: alignment.sam
        samtools stats alignment.sam > samtools_stats.txt

        if [ -f "alignment.sam" ] && [ -f "samtools_stats.txt" ]
        then
            echo pass > .status
        else
            echo fail > .status
        fi

        echo -n "" > .report.json
        echo -n "" > .versions
    } || {
        echo fail > .status
    }
    """
}


process insert_size_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/assembly/insert_size_{{ pid }}/"

    input:
    set sample_id, file(sam_stats) from IN_insert_size_{{ pid }}
    val plot from IN_plot{{ pid }}

    output:
    file ("*insert_size_report.tab")
    file ("*insert_size_distribution.html") optional true
    {% with task_name="insert_size" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "insert_size.py"

}