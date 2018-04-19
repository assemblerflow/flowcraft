
process pilon_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    echo false
    publishDir 'results/assembly/pilon_{{ pid }}/', mode: 'copy', pattern: "*.fasta"

    input:
    set sample_id, file(assembly), file(bam_file), file(bam_index) from {{ input_channel }}

    output:
    set sample_id, '*_polished.fasta' into {{ output_channel }}, pilon_report_{{ pid }}
    {% with task_name="pilon" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        pilon_mem=${String.valueOf(task.memory).substring(0, String.valueOf(task.memory).length() - 1).replaceAll("\\s", "")}
        java -jar -Xms256m -Xmx\${pilon_mem} /NGStools/pilon-1.22.jar --genome $assembly --frags $bam_file --output ${assembly.name.replaceFirst(~/\.[^\.]+$/, '')}_polished --changes --threads $task.cpus >> .command.log 2>&1
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}

process pilon_report_{{ pid }} {

    {% with overwrite="false" %}
    {% include "report_post.txt" ignore missing %}
    {% endwith %}

    tag { sample_id }

    input:
    set sample_id, file(assembly), file(coverage_bp) from pilon_report_{{ pid }}.join(SIDE_BpCoverage_{{ pid }})

    output:
    file "*_assembly_report.csv" into pilon_report_out_{{ pid }}
    {% with task_name="pilon_report" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "assembly_report.py"

}


process compile_pilon_report_{{ pid }} {

    publishDir "reports/assembly/pilon_{{ pid }}/", mode: 'copy'

    input:
    file(report) from pilon_report_out_{{ pid }}.collect()

    output:
    file "pilon_assembly_report.csv"

    """
    echo Sample,Number of contigs,Average contig size,N50,Total assembly length,GC content,Missing data > pilon_assembly_report.csv
    cat $report >> pilon_assembly_report.csv
    """
}

{{ forks }}

