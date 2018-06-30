file(params.referenceFileO{{ param_id }}) ? params.referenceFileO{{ param_id }} : exit(1, "'referenceFileO{{ param_id }}' parameter missing")
IN_refO_{{ pid }} = Channel.fromPath(params.referenceFileO{{ param_id }})
    .map{ it -> it.exists() ? it : exit(1, "referenceFileO file was not found: '${params.referenceFileO{{ param_id }}}'")}

file(params.referenceFileH{{ param_id }}) ? params.referenceFileH{{ param_id }} : exit(1, "'referenceFileH{{ param_id }}' parameter missing")
IN_refH_{{ pid }} = Channel.fromPath(params.referenceFileH{{ param_id }})
    .map{ it -> it.exists() ? it : exit(1, "referenceFileH file was not found: '${params.referenceFileH{{ param_id }}}'")}

process seq_typing_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    errorStrategy "ignore"
    publishDir "results/seqtyping/${sample_id}/"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    file refO from IN_refO_{{ pid }}
    file refH from IN_refH_{{ pid }}

    output:
    file "seq_typing*"
    {% with task_name="seq_typing" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # Prevents read-only issues
        mkdir rematch_temp
        cp -r /NGStools/ReMatCh rematch_temp
        export PATH="\$(pwd)/rematch_temp/ReMatCh:\$PATH"

        seq_typing.py -f ${fastq_pair[0]} ${fastq_pair[1]} -r \$(pwd)/$refO \$(pwd)/$refH -o ./ -j $task.cpus --extraSeq 0 --mapRefTogether --minGeneCoverage 60
        json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'seqtyping','value':'\$(cat seq_typing.report.txt)','table':'typing'}]}]}"
        echo \$json_str > .report.json

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

