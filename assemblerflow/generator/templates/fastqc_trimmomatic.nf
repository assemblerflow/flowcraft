
process fastqc {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }

    input:
    set fastq_id, file(fastq_pair) from {{ input_channel }}
    val ad from Channel.value('None')

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into MAIN_fastqc_out_{{ pid }}
    set fastq_id, val("fastqc_{{ pid }}"), file(".status"), file(".warning"), file(".fail") into STATUS_fastqc_{{ pid }}
    file ".report.json"

    when:
    params.stopAt != "fastqc"

    script:
    template "fastqc.py"
}

/** FASTQC_REPORT - MAIN
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report {

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    tag { fastq_id + " getStats" }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc/run_1/', pattern: '*summary.txt', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_{{ pid }}
    val opts from Channel.value("--ignore-tests")

    output:
    set fastq_id, file(fastq_pair), 'optimal_trim', ".status" into MAIN_fastqc_trim
    file '*_trim_report' into LOG_trim_{{ pid }}
    set fastq_id, val("fastqc_report_{{ pid }}"), file(".status"), file(".warning"), file(".fail") into STATUS_report_{{ pid }}
    file "*_status_report" into LOG_fastqc_report_{{ pid }}
    file "${fastq_id}_*_summary.txt" optional true
    file ".report.json"

    script:
    template "fastqc_report.py"

}

MAIN_fastqc_trim_{{ pid }} = Channel.create()
MAIN_fastqc_trim
        .filter{ it[3].text == "pass" }
        .map{ [it[0], it[1], file(it[2]).text] }
        .into(MAIN_fastqc_trim_{{ pid }})


/** TRIM_REPORT - PLUG-IN
This will collect the optimal trim points assessed by the fastqc_report
process and write the results of all samples in a single csv file
*/
process trim_report {

    publishDir 'reports/fastqc/', mode: 'copy'

    input:
    file trim from LOG_trim_{{ pid }}.collect()

    output:
    file "FastQC_trim_report.csv"

    """
    echo Sample,Trim begin, Trim end >> FastQC_trim_report.csv
    cat $trim >> FastQC_trim_report.csv
    """
}


process compile_fastqc_status {

    publishDir 'reports/fastqc/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_{{ pid }}.collect()

    output:
    file 'FastQC_1run_report.csv'

    """
    echo Sample, Failed? >> FastQC_1run_report.csv
    cat $rep >> FastQC_1run_report.csv
    """

}


/** TRIMMOMATIC - MAIN
This process will execute trimmomatic. Currently, the main channel requires
information on the trim_range and phred score.
*/
process trimmomatic {

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    tag { fastq_id + " getStats" }

    input:
    set fastq_id, file(fastq_pair), trim_range, phred from MAIN_fastqc_trim_{{ pid }}.join(SIDE_phred_{{ pid }})
    val opts from IN_trimmomatic_opts

    output:
    set fastq_id, "${fastq_id}_*P*" optional true into {{ output_channel }}
    set fastq_id, val("trimmomatic_{{ pid }}"), file(".status"), file(".warning"), file(".fail") into STATUS_trim_{{ pid }}
    file 'trimmomatic_report.csv'
    file ".report.json"

    when:
    params.stopAt != "trimmomatic"

    script:
    template "trimmomatic.py"

}

{{ forks }}

