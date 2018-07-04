// Check sliding window parameter
if ( params.trimSlidingWindow{{ param_id }}.toString().split(":").size() != 2 ){
    exit 1, "'trimSlidingWindow{{ param_id }}' parameter must contain two values separated by a ':'. Provided value: '${params.trimSlidingWindow{{ param_id }}}'"
}
if ( !params.trimLeading{{ param_id }}.toString().isNumber() ){
    exit 1, "'trimLeading{{ param_id }}' parameter must be a number. Provide value: '${params.trimLeading{{ param_id }}}'"
}
if ( !params.trimTrailing{{ param_id }}.toString().isNumber() ){
    exit 1, "'trimTrailing{{ param_id }}' parameter must be a number. Provide value: '${params.trimTrailing{{ param_id }}}'"
}
if ( !params.trimMinLength{{ param_id }}.toString().isNumber() ){
    exit 1, "'trimMinLength{{ param_id }}' parameter must be a number. Provide value: '${params.trimMinLength{{ param_id }}}'"
}

IN_trimmomatic_opts_{{ pid }} = Channel.value([params.trimSlidingWindow{{ param_id }},params.trimLeading{{ param_id }},params.trimTrailing{{ param_id }},params.trimMinLength{{ param_id }}])
IN_adapters_{{ pid }} = Channel.value(params.adapters{{ param_id }})

clear = params.clearAtCheckpoint ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process fastqc_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "reports/fastqc_{{ pid }}/", pattern: "*.html"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val ad from Channel.value('None')

    output:
    set sample_id, file(fastq_pair), file('pair_1*'), file('pair_2*') into MAIN_fastqc_out_{{ pid }}
    file "*html"
    {% with task_name="fastqc" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "fastqc.py"
}

/** FASTQC_REPORT - MAIN
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report_{{ pid }} {

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc_{{ pid }}/run_1/', pattern: '*summary.txt', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_{{ pid }}
    val opts from Channel.value("--ignore-tests")

    output:
    set sample_id, file(fastq_pair), 'optimal_trim', ".status" into _MAIN_fastqc_trim_{{ pid }}
    file '*_trim_report' into LOG_trim_{{ pid }}
    file "*_status_report" into LOG_fastqc_report_{{ pid }}
    file "${sample_id}_*_summary.txt" optional true
    {% with task_name="fastqc_report" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "fastqc_report.py"

}

MAIN_fastqc_trim_{{ pid }} = Channel.create()
_MAIN_fastqc_trim_{{ pid }}
        .filter{ it[3].text == "pass" }
        .map{ [it[0], it[1], file(it[2]).text] }
        .into(MAIN_fastqc_trim_{{ pid }})


/** TRIM_REPORT - PLUG-IN
This will collect the optimal trim points assessed by the fastqc_report
process and write the results of all samples in a single csv file
*/
process trim_report_{{ pid }} {

    publishDir 'reports/fastqc_{{ pid }}/', mode: 'copy'

    input:
    file trim from LOG_trim_{{ pid }}.collect()

    output:
    file "FastQC_trim_report.csv"

    """
    echo Sample,Trim begin, Trim end >> FastQC_trim_report.csv
    cat $trim >> FastQC_trim_report.csv
    """
}


process compile_fastqc_status_{{ pid }} {

    publishDir 'reports/fastqc_{{ pid }}/', mode: 'copy'

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
process trimmomatic_{{ pid }} {

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    tag { sample_id }
    publishDir "results/trimmomatic_{{ pid }}", pattern: "*.gz"

    input:
    set sample_id, file(fastq_pair), trim_range, phred from MAIN_fastqc_trim_{{ pid }}.join(SIDE_phred_{{ pid }})
    val opts from IN_trimmomatic_opts_{{ pid }}
    val ad from IN_adapters_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, "${sample_id}_*trim.fastq.gz" into {{ output_channel }}
    file 'trimmomatic_report.csv'
    {% with task_name="trimmomatic" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "trimmomatic.py"

}

{{ forks }}

