
process momps_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})

    output:
    file("*_st.tsv") into momps_st_{{ pid }}
    file("*_profile.tsv") into momps_profile_{{ pid }}
    {% with task_name="momps" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # Stage in momps source files. This cannot be a symlink because the files
        # need to be writable.
        cp -r /NGStools/mompS/* .
        momps.pl -r ${fastq[0]} -f ${fastq[1]} -a $assembly -o res -p $sample_id -t ${task.cpus}
        # Get the ST for the sample
        if [ -f "res/${sample_id}.MLST_res.txt" ]
        then
            st=\$(grep -oP "ST = \\K\\w+" res/*.MLST_res.txt)
            # If the ST cannot be determined, set string to ND
            if [ -z \$st ]
            then
                st="ND"
            fi
            echo $sample_id\t\${st}> ${sample_id}_st.tsv
            # Add ST information to report JSON
            json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'mompS','value':'\$st','table':'typing'}]}]}"
            echo \$json_str > .report.json
            # Get the profile for the sample
            echo $sample_id\t\$(awk "NR == 7" res/*.MLST_res.txt) > ${sample_id}_profile.tsv
            rm -r res
        else
            echo fail > .status
            rm -r res
        fi
    } || {
        echo fail > .status
        # Remove results directory
        rm -r res
    }
    """

}


process momps_report_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}
    publishDir "results/typing/momps_{{ pid }}/", pattern: "*.tsv"

    input:
    file(st_file) from momps_st_{{ pid }}.collect()
    file(profile_file) from momps_profile_{{ pid }}.collect()

    output:
    file "*.tsv"


    script:
    """
    cat $st_file >> momps_st.tsv
    cat $profile_file >> momps_profile.tsv
    """

}

{{ forks }}

