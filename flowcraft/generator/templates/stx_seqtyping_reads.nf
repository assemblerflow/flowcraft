stx2covered = ""
stx2identity = ""

if ( params.stx2covered{{ param_id }} ){
  if ( ! params.stx2covered{{ param_id }}.toString().isNumber() ){
    exit 1, "--'stx2covered{{ param_id }}' parameter must be a number. Provided value: '${params.stx2covered{{ param_id }}}'"
  }
  stx2covered = "--stx2covered ${params.stx2covered{{ param_id }}}"
}

if ( params.stx2identity{{ param_id }} ){
  if ( ! params.stx2identity{{ param_id }}.toString().isNumber() ){
    exit 1, "--'stx2identity{{ param_id }}' parameter must be a number. Provided value: '${params.stx2identity{{ param_id }}}'"
  }
  stx2identity = "--stx2identity ${params.stx2identity{{ param_id }}}"
}

process stx_seqtyping_reads_{{ pid }} {
    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    errorStrategy { task.exitStatus == 120 ? 'ignore' : 'ignore' }
    publishDir path: "results/typing/stx_seqtyping_reads/${sample_id}/", mode: 'symlink', overwrite: true, pattern: 'seq_typing.ecoli_stx_subtyping.*'

    input:
    set sample_id, file(fastq) from {{ input_channel }}

    output:
    file "seq_typing.ecoli_stx_subtyping.*"
    {% with task_name="stx_seqtyping_reads" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    exit_code=0

    version_str="[{'program':'ecoli_stx_subtyping.py','version':'\$(ecoli_stx_subtyping.py --version | cut -d \' \' -f 2)'}]"
    echo \$version_str > .versions

    status='error'
    report_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'stx_type_seqtyping_reads','value':'NA','table':'typing'}]}]}"

    {
      ecoli_stx_subtyping.py reads -f $fastq --org stx subtyping -o ./ -j $task.cpus $stx2covered $stx2identity
    } || {
      exit_code=\$?
    }

    if [ \$exit_code -eq 0 ]; then
      status='pass'

      stx_type=\$(cat seq_typing.ecoli_stx_subtyping.txt)

      report_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'stx_type_seqtyping_reads','value':'\$stx_type','table':'typing'}]}]}"
    fi

    echo \$status > .status
    echo \$report_str > .report.json

    exit \$exit_code
    """
}

{{ forks }}
