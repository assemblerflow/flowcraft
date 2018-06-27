IN_genome_size_{{ pid }} = Channel.value(params.genomeSize{{ param_id }})
    .map{it -> it.toString().isNumber() ? it : exit(1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize_{{ param_id }}}'")}

IN_min_coverage_{{ pid }} = Channel.value(params.minCoverage{{ param_id }})
    .map{it -> it.toString().isNumber() ? it : exit(1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage_{{ param_id }}}'")}

process integrity_coverage_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val gsize from IN_genome_size_{{ pid }}
    val cov from IN_min_coverage_{{ pid }}
    // This channel is for the custom options of the integrity_coverage.py
    // script. See the script's documentation for more information.
    val opts from Channel.value('')

    output:
    set sample_id,
        file(fastq_pair),
        file('*_encoding'),
        file('*_phred'),
        file('*_coverage'),
        file('*_max_len') into MAIN_integrity_{{ pid }}
    file('*_report') optional true into LOG_report_coverage1_{{ pid }}
    {% with task_name="integrity_coverage" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
LOG_corrupted_{{ pid }} = Channel.create()
MAIN_PreCoverageCheck_{{ pid }} = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
MAIN_integrity_{{ pid }}.choice(LOG_corrupted_{{ pid }}, MAIN_PreCoverageCheck_{{ pid }}) {
    a -> a[2].text == "corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
{{ output_channel }} = Channel.create()
SIDE_phred_{{ pid }} = Channel.create()
SIDE_max_len_{{ pid }} = Channel.create()

MAIN_PreCoverageCheck_{{ pid }}
// Low coverage samples have the 4th value of the Channel with 'fail'
    .filter{ it[4].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'SIDE_phred'
    .separate({{ output_channel }}, SIDE_phred_{{ pid }}, SIDE_max_len_{{ pid }}){
        a -> [ [a[0], a[1]], [a[0], a[3].text], [a[0], a[5].text]  ]
    }

/** REPORT_COVERAGE - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_initial.csv'
*/
process report_coverage_{{ pid }} {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage_{{ pid }}/'

    input:
    file(report) from LOG_report_coverage1_{{ pid }}.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_initial.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_initial.csv
    cat $report >> estimated_coverage_initial.csv
    """
}

/** REPORT_CORRUPT - PLUG-IN
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt_{{ pid }} {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/corrupted_{{ pid }}/'

    input:
    val sample_id from LOG_corrupted_{{ pid }}.collect{it[0]}

    output:
    file 'corrupted_samples.txt'

    """
    echo ${sample_id.join(",")} | tr "," "\n" >> corrupted_samples.txt
    """

}

{{ forks }}

