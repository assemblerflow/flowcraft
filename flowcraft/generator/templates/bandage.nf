if (params.reference{{param_id}} != null) {
    process bandage_{{pid}} {
        {% include "post.txt" ignore missing %}

        tag { sample_id }
        publishDir "reports/assembly/bandage_{{pid}}/$sample_id"

        input:
        set sample_id, file(assembly) from {{input_channel}}
        file reference from Channel.fromPath(params.reference{{param_id}})

        output:
        file "*.png *.svg"
        {% with task_name="bandage" %}
        {%- include "compiler_channels.txt" ignore missing -%}
        {% endwith %}

        script:
        """
        time Bandage image $assembly ${assembly}.png >>.command.log 2>&1
        time Bandage image $assembly ${assembly}.svg >>.command.log 2>&1
        time Bandage image $assembly ${assembly}.ref.png --query $reference >>.command.log 2>&1
        time Bandage image $assembly ${assembly}.ref.svg --query $reference >>.command.log 2>&1
        """
    }
} else {
    process bandage_{{pid}} {
        {% include "post.txt" ignore missing %}

        tag { sample_id }
        publishDir "reports/assembly/bandage_{{pid}}/$sample_id"

        input:
        set sample_id, file(assembly) from {{input_channel}}

        output:
        file "*.png *.svg"
        {% with task_name="bandage" %}
        {%- include "compiler_channels.txt" ignore missing -%}
        {% endwith %}

        script:
        """
        time Bandage image $assembly ${assembly}.png >>.command.log 2>&1
        time Bandage image $assembly ${assembly}.svg >>.command.log 2>&1
        """
    }
}
