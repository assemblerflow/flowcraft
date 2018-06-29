IN_adapters_{{ pid }} = Channel.value(params.adapters{{ param_id }})

process fastqc2_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "reports/fastqc_{{ pid }}/", pattern: "*.html"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val ad from IN_adapters_{{ pid }}

    output:
    set sample_id, file(fastq_pair), file('pair_1*'), file('pair_2*') into MAIN_fastqc_out_{{ pid }}
    file "*html"
    {% with task_name="fastqc2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "fastqc.py"
}


process fastqc2_report_{{ pid }} {

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc_{{ pid }}/run_2/', pattern: '*summary.txt', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_{{ pid }}
    val opts from Channel.value("")

    output:
    set sample_id, file(fastq_pair), '.status' into MAIN_fastqc_report_{{ pid }}
    file "*_status_report" into LOG_fastqc_report_{{ pid }}
    file "${sample_id}_*_summary.txt" optional true
    {% with task_name="fastqc2_report" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "fastqc_report.py"

}


process compile_fastqc_status2_{{ pid }} {

    publishDir 'reports/fastqc_{{ pid }}/', mode: 'copy'

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

