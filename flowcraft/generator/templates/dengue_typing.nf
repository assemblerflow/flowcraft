process dengue_typing_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    errorStrategy { task.exitStatus == 120 ? 'ignore' : 'ignore' }
    publishDir "results/dengue_typing/${sample_id}/"


    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}

    output:
    file "seq_typing*"
    set sample_id, file(assembly) into {{ output_channel}}
    file("*.fa") optional true into _ref_seqTyping_{{ pid }}
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

        echo \$(cat ${assembly})

        if [ \$(cat ${assembly}) == "fail" ];
        then
            seq_typing.py reads --org Dengue Virus  -o ./ -j $task.cpus -f ${fastq[0]} ${fastq[1]}

            json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'seqtyping','value':'\$(cat seq_typing.report.txt)','table':'typing'}]}],'column':'typing'}]}"

            rm ${assembly} && exit 120

        else
            seq_typing.py assembly --org Dengue Virus -f ${assembly} -o ./ -j $task.cpus -t nucl

            if [ ${params.reference{{ param_id}} } == "true" ];
            then
                awk 'NR == 2 { print \$4 }' seq_typing.report_types.tab > reference
                parse_fasta.py -t \$(cat reference)  -f /NGStools/seq_typing/reference_sequences/dengue_virus/1_GenotypesDENV_14-05-18.fasta
                json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'seqtyping','value':'\$(cat seq_typing.report.txt)','table':'typing'}]}],'metadata':[{'sample':'${sample_id}','treeData':'\$(cat seq_typing.report.txt)','column':'typing'},{'sample':'\$(cat header.txt)','treeData':'\$(cat seq_typing.report.txt)','column':'typing'}]}"
            else
                json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'seqtyping','value':'\$(cat seq_typing.report.txt)','table':'typing'}]}],'metadata':[{'sample':'${sample_id}','treeData':'\$(cat seq_typing.report.txt)','column':'typing'}]}"
            fi

        fi


        # Add information to dotfiles
        echo \$json_str > .report.json
        version_str="[{'program':'seq_typing.py','version':'0.1'}]"
        echo \$version_str > .versions
        rm -r rematch_temp
        if [ -s seq_typing.report.txt ];
        then
            echo pass > .status
        else
            echo fail > .status
        fi
    } || {
        echo fail > .status
        json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'seqtyping','value':'NA','table':'typing'}]}]}"
        echo \$json_str > .report.json
    }
    """

}

{{ forks }}

