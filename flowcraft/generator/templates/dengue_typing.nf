file(params.BD_sequence_file{{ param_id }}) ? params.BD_sequence_file{{ param_id }} : exit(1, "'BD_sequence_file{{ param_id }}' parameter missing")

process dengue_typing_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    errorStrategy "ignore"
    publishDir "results/dengue_typing/${sample_id}/"

    input:
    set sample_id, file(assembly) from {{ input_channel }}

    output:
    file "seq_typing*"
    {% with task_name="dengue_typing" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # Prevents read-only issues
        mkdir rematch_temp
        cp -r /NGStools/ReMatCh rematch_temp
        export PATH="\$(pwd)/rematch_temp/ReMatCh:\$PATH"

        seq_typing.py assembly -f ${assembly} -b ${ params.BD_sequence_file{{ param_id }} } -o ./ -j $task.cpus -t nucl

        rm -r rematch_temp

        if [ -s seq_typing.report.txt ];
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

