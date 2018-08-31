
process sistr_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/typing/sistr_{{ pid }}', pattern: ".tab", mode: "copy"

    input:
    set sample_id, file(assembly) from {{ input_channel }}

    output:
    {% with task_name="sistr" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        sistr --qc -vv -t $task.cpus -f tab -o ${sample_id}_sistr.tab ${assembly}
        json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'sistr','value':'\$(awk \"FNR == 2\" *.tab | cut -f14)','table':'typing'}]}]}"
        echo \$json_str > .report.json
        sistr_version=\$(sistr --version | cut -d" " -f2)
        version_str="[{'program':'sistr','version':'\$sistr_version'}]"
        echo \$version_str > .versions

        if [ -s ${sample_id}_sistr.tab ];
        then
            echo pass > .status
        else
            echo fail > .status
        fi

    } || {
        echo fail > .status
    }
    """
}

{{ forks }}
