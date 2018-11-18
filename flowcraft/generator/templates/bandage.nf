// True when a GFA secondary channel is connected to this component.
has_gfa1_{{pid}} = binding.hasVariable('gfa1_{{pid}}')

process bandage_{{pid}} {
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "reports/assembly/bandage_{{pid}}/$sample_id"

    input:
    set sample_id, file(fasta) from {{input_channel}}
    file gfa1 from has_gfa1_{{pid}} ? gfa1_{{pid}} : Channel.value("NA")
    file reference from params.reference{{param_id}} ?
        Channel.fromPath(params.reference{{param_id}}) :
        Channel.value("NA")

    output:
    file "*.png"
    file "*.svg"
    {% with task_name="bandage" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    // Use the GFA assembly when available and FASTA otherwise.
    assembly = has_gfa1_{{pid}} ? gfa1 : fasta
    command =
        """
        time Bandage image $assembly ${assembly}.png >>.command.log 2>&1
        time Bandage image $assembly ${assembly}.svg >>.command.log 2>&1
        """
    if (params.reference{{param_id}})
        command +=
            """
            time Bandage image $assembly ${assembly}.ref.png --query $reference >>.command.log 2>&1
            time Bandage image $assembly ${assembly}.ref.svg --query $reference >>.command.log 2>&1
            """
    command
}
