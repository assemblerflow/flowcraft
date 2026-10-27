// True when a raw_long_reads secondary channel is connected to this component.
has_raw_long_reads_{{pid}} = binding.hasVariable('raw_long_reads_{{pid}}')

process porechop_{{pid}} {
    {% include "post.txt" ignore missing %}

    publishDir "results/porechop_{{pid}}", pattern: "*.fastq.gz"
    publishDir "reports/porechop_{{pid}}", pattern: "*.log"

    tag { sample_id }

    input:
    set sample_id, file(fastq_pair) from {{input_channel}}
    file raw_long_reads from has_raw_long_reads_{{pid}} ? raw_long_reads_{{pid}} :
        Channel.fromPath(params.long_reads{{param_id}})

    output:
    set sample_id, file(fastq_pair) into {{output_channel}}
    file "${sample_id}.porechop.fastq.gz" into long_reads_{{pid}}
    {% with task_name="porechop" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    "time porechop -t $task.cpus --format fastq.gz -i ${raw_long_reads} -o ${sample_id}.porechop.fastq.gz >${sample_id}.log"
}

{{ forks }}
