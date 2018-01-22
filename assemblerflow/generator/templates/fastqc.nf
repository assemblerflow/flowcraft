
process fastqc2 {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }

    input:
    set fastq_id, file(fastq_pair) from {{ input_channel }}
    val ad from IN_adapters

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into MAIN_fastqc_out_{{ pid }}
    set fastq_id, val("fastqc2_{{ pid }}"), file(".status"), file(".warning"), file(".fail") into STATUS_fastqc_{{ pid }}
    file ".report.json"

    when:
    params.stopAt != "fastqc2"

    script:
    template "fastqc.py"
}


process fastqc2_report {

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    tag { fastq_id + " getStats" }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc/run_2/', pattern: '*summary.txt', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_{{ pid }}
    val opts from Channel.value("")

    output:
    set fastq_id, file(fastq_pair), '.status' into MAIN_fastqc_report_{{ pid }}
    file "*_status_report" into LOG_fastqc_report_{{ pid }}
    set fastq_id, val("fastqc2_report_{{ pid }}"), file(".status"), file(".warning"), file(".fail") into STATUS_report_{{ pid }}
    file "${fastq_id}_*_summary.txt" optional true
    file ".report.json"

    script:
    template "fastqc_report.py"

}


process compile_fastqc_status2 {

    publishDir 'reports/fastqc/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_{{ pid }}.collect()

    output:
    file 'FastQC_2run_report.csv'

    """
    echo Sample, Failed? >> FastQC_2run_report.csv
    cat $rep >> FastQC_2run_report.csv
    """

}

{{ output_channel }} = Channel.create()

MAIN_fastqc_report_{{ pid }}
        .filter{ it[2].text == "pass" }
        .map{ [it[0], it[1]] }
        .into({{ output_channel }})

{{ forks }}

