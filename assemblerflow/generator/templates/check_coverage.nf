
process integrity_coverage_2 {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    cpus 1

	input:
	set fastq_id, file(fastq_pair) from {{ input_channel }}
	val gsize from IN_genome_size
	val cov from IN_min_coverage
	// Use -e option for skipping encoding guess
	val opts from Channel.value('-e')

	output:
	set fastq_id,
	    file(fastq_pair),
	    file('*_coverage'),
	    file('*_max_len') optional true into MAIN_integrity_{{ pid }}
	file('*_report') into LOG_report_coverage_{{ pid }}
	set fastq_id, val("integrity_coverage2_{{ pid }}"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}
	file ".report.json"

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


process report_coverage_2 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage/'

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

