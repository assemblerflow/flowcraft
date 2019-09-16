IN_min_contig_{{ pid }} = Channel.value(params.minContig{{ param_id }})
IN_min_align_score_{{ pid }} = Channel.value(params.minAlignScore{{ param_id }})

clear = params.clearInput{{ param_id }} ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process vamb_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    //publishDir "results/assembly/binning/vamb_{{ pid }}/${sample_id}/"

    input:
    set sample_id, file(assembly), file(bam_file), file(bam_index) from {{ input_channel }}
    val length_threshold from IN_min_contig_{{ pid }}
    val min_score from IN_min_align_score_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:

    {% with task_name="vamb"%}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # run METABAT2
        run.py results/ ${assembly} ${bam_file} -m ${length_threshold} -s ${min_score}

        # In case no sequences are binned
        if [ -z "\$(ls -A *metabat-bins*/)" ]; then
            echo "false" > false_bin.fa
            mv false_bin.fa *metabat-bins*/
            echo "false" > bin_status.txt;
        else
            echo "true" > bin_status.txt
        fi

    } || {
        echo fail > .status
    }
    """
}


{{ forks }}