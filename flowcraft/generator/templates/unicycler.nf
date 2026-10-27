// True when a long_reads secondary channel is connected to this component.
has_long_reads_{{pid}} = binding.hasVariable('long_reads_{{pid}}')

process unicycler_{{pid}} {
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/unicycler_{{pid}}/', pattern: 'assembly.fasta'
    publishDir 'results/assembly/unicycler_{{pid}}/', pattern: 'assembly.gfa'

    input:
    set sample_id, file(fastq_pair) from {{input_channel}}
    file long_reads from has_long_reads_{{pid}} ? long_reads_{{pid}} :
        params.long_reads{{param_id}} ? Channel.fromPath(params.long_reads{{param_id}}) :
        Channel.value("NA")

    output:
    set sample_id, file('assembly.fasta') into {{output_channel}}
    file "assembly.gfa" into gfa1_{{pid}}
    {% with task_name="unicycler" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    command = "unicycler -t $task.cpus -o . --no_correct --no_pilon -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}"
    if (params.long_reads{{param_id}})
        command += " -l ${long_reads}"
    command
}

{{forks}}
