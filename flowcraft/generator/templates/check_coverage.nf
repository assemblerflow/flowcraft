IN_genome_size_{{ pid }} = Channel.value(params.genomeSize{{ param_id }})
    .map{it -> it.toString().isNumber() ? it : exit (1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize{{ param_id }}}'")}
IN_min_coverage_{{ pid }} = Channel.value(params.minCoverage{{ param_id }})
    .map{it -> it.toString().isNumber() ? it : exit (1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage{{ param_id }}}'")}

process integrity_coverage2_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    cpus 1

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val gsize from IN_genome_size_{{ pid }}
    val cov from IN_min_coverage_{{ pid }}
    // Use -e option for skipping encoding guess
    val opts from Channel.value('-e')

    output:
    set sample_id,
        file(fastq_pair),
        file('*_coverage'),
        file('*_max_len') optional true into MAIN_integrity_{{ pid }}
    file('*_report') into LOG_report_coverage_{{ pid }}
    {% with task_name="check_coverage" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "integrity_coverage.py"
}

{{ output_channel }} = Channel.create()
SIDE_max_len_{{ pid }} = Channel.create()

MAIN_integrity_{{ pid }}
    .filter{ it[2].text != "fail" }
    .separate({{ output_channel }}, SIDE_max_len_{{ pid }}){
        a -> [ [a[0], a[1]], [a[0], a[3].text]]
    }


process report_coverage2_{{ pid }} {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage_{{ pid }}/'

    input:
    file(report) from LOG_report_coverage_{{ pid }}.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_second.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_second.csv
    cat $report >> estimated_coverage_second.csv
    """
}

{{ forks }}

