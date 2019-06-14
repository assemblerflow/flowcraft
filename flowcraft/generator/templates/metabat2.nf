IN_contig_percentage_{{ pid }} = Channel.value(params.maxPercentage{{ param_id }})
IN_length_threshold_{{ pid }} = Channel.value(params.minContig{{ param_id }})

clear = params.clearInput{{ param_id }} ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process metabat2_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    //publishDir "results/assembly/binning/metabat2_{{ pid }}/${sample_id}/"

    input:
    set sample_id, file(assembly), file(bam_file), file(bam_index) from {{ input_channel }}
    val contig_percentage from IN_contig_percentage_{{ pid }}
    val length_threshold from IN_length_threshold_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, file(assembly), file('*metabat-bins*/*.fa'), file ('bin_status.txt') into binCh_{{ pid }}
    set sample_id, file('*metabat-bins*/*.fa') into intoReport_{{ pid }}
    {% with task_name="metabat2"%}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # prevent indexing errors
        samtools sort ${bam_file} sorted
        samtools index sorted.bam

        # run METABAT2
        runMetaBat.sh -m ${length_threshold} --unbinned --maxP ${contig_percentage} ${assembly} sorted.bam

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

process report_metabat2_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(bins) from intoReport_{{ pid }}

    output:
    {% with task_name="report_metabat2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_metabat.py"
}

// If maxbin fails to obtain bins for a sample, the workflow continues with the original assembly
{{ output_channel }} = Channel.create()

OUT_binned = Channel.create()
OUT_unbinned = Channel.create()

failedBinning = Channel.create()
successfulBinning = Channel.create()

binCh_{{ pid }}.choice(failedBinning, successfulBinning){ it -> it[3].text == "false\n" ? 0 : 1 }

failedBinning.map{ it -> [it[0], it[1]] }.into(OUT_unbinned)

successfulBinning.map{ it -> [it[2].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it[2]]}
    .transpose()
    .map{it -> [it[1].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'),it[1]]}
    .into(OUT_binned)

OUT_binned.mix(OUT_unbinned).set{ {{ output_channel }} }

{{ forks }}