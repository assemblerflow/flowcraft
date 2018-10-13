if (params.reference{{param_id}} == null && params.genomeSizeBp{{param_id}} == null)
    exit 1, "Specify at least one of reference or genomeSizeBp"
if (params.reference{{param_id}} != null && params.genomeSizeBp{{param_id}} != null)
    exit 1, "Specify only one of reference or genomeSizeBp"

if (params.reference{{param_id}} != null) {
    process quast_{{pid}} {
        {% include "post.txt" ignore missing %}

        tag { sample_id }
        publishDir "results/assembly/quast_{{pid}}/$sample_id", pattern: "*.tsv"
        publishDir "reports/assembly/quast_{{pid}}/$sample_id"

        input:
        set sample_id, file(assembly) from {{input_channel}}
        file reference from Channel.fromPath(params.reference{{param_id}})

        output:
        file "*"
        {% with task_name="quast" %}
        {%- include "compiler_channels.txt" ignore missing -%}
        {% endwith %}

        script:
        "/usr/bin/time -v quast -o . -r $reference -s $assembly -l $sample_id -t $task.cpus >> .command.log 2>&1"
    }
} else if (params.genomeSizeBp{{param_id}} != null) {
    process quast_{{pid}} {
        {% include "post.txt" ignore missing %}

        tag { sample_id }
        publishDir "results/assembly/quast_{{pid}}/$sample_id", pattern: "*.tsv"
        publishDir "reports/assembly/quast_{{pid}}/$sample_id"

        input:
        set sample_id, file(assembly) from {{input_channel}}
        val genomeSizeBp from Channel.value(params.genomeSizeBp{{param_id}})

        output:
        file "*"
        {% with task_name="quast" %}
        {%- include "compiler_channels.txt" ignore missing -%}
        {% endwith %}

        script:
        "/usr/bin/time -v quast -o . --est-ref-size=$genomeSizeBp -s $assembly -l $sample_id -t $task.cpus >> .command.log 2>&1"
    }
}
