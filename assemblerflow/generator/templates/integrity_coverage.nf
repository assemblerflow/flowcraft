
process integrity_coverage {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    // This process can only use a single CPU
    cpus 1

	input:
	set fastq_id, file(fastq_pair) from {{ input_channel }}
	val gsize from IN_genome_size
	val cov from IN_min_coverage
	// This channel is for the custom options of the integrity_coverage.py
	// script. See the script's documentation for more information.
	val opts from Channel.value('')

	output:
	set fastq_id,
	    file(fastq_pair),
	    file('*_encoding'),
	    file('*_phred'),
	    file('*_coverage'),
	    file('*_max_len') optional true into MAIN_integrity
	file('*_report') optional true into LOG_report_coverage1
	set fastq_id, val("integrity_coverage_{{ pid }}"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}
	file ".report.json"

	script:
	template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
LOG_corrupted = Channel.create()
MAIN_PreCoverageCheck = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
MAIN_integrity.choice(LOG_corrupted, MAIN_PreCoverageCheck) {
    a -> a[2].text == "corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
{{ output_channel }} = Channel.create()
SIDE_phred_{{ pid }} = Channel.create()
SIDE_max_len_{{ pid }} = Channel.create()

MAIN_PreCoverageCheck
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
process report_coverage {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage/'

    input:
    file(report) from LOG_report_coverage1.filter{ it.text != "corrupt" }.collect()

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
process report_corrupt {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/corrupted/'

    input:
    val fastq_id from LOG_corrupted.collect{it[0]}

    output:
    file 'corrupted_samples.txt'

    """
    echo ${fastq_id.join(",")} | tr "," "\n" >> corrupted_samples.txt
    """

}

{{ forks }}

