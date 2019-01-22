IN_min_contig_lenght_{{ pid }} = Channel.value(params.min_contig_lenght{{ param_id }})
IN_max_iteration_{{ pid }} = Channel.value(params.max_iteration{{ param_id }})
IN_prob_threshold_{{ pid }} = Channel.value(params.prob_threshold{{ param_id }})

clear = params.clearInput{{ param_id }} ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process maxbin2_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/maxbin2_{{ pid }}/${sample_id}/"

    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})
    val minContigLenght from IN_min_contig_lenght_{{ pid }}
    val maxIterations from IN_max_iteration_{{ pid }}
    val probThreshold from IN_prob_threshold_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, file(assembly), file ('*_maxbin.*.fasta'), file ('bin_status.txt') into binCh_{{ pid }}
    file '*_maxbin.{abundance,log,summary}'
    set sample_id, file("*_maxbin.summary") into intoReport_{{ pid }}

    {% with task_name="maxbin2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        run_MaxBin.pl -contig ${assembly} -out ${sample_id}_maxbin -reads ${fastq[0]} -reads2 ${fastq[1]} \
        -thread $task.cpus -min_contig_length ${minContigLenght} -max_iteration ${maxIterations} \
        -prob_threshold ${probThreshold}

        echo pass > .status

        #in case maxbin fails to bin sequences for a sample:
        if ls *_maxbin.*.fasta 1> /dev/null 2>&1; then echo "true" > bin_status.txt; else echo "false" \
        > false_maxbin.0.fasta; echo "false" > bin_status.txt; fi


        if [ "$clear" = "true" ];
        then
            work_regex=".*/work/.{2}/.{30}/.*"
            file_source1=\$(readlink -f \$(pwd)/${fastq[0]})
            file_source2=\$(readlink -f \$(pwd)/${fastq[1]})
            assembly_file=\$(readlink -f \$(pwd)/${assembly})
            if [[ "\$file_source1" =~ \$work_regex ]]; then
                rm \$file_source1 \$file_source2 \$assembly_file
            fi
        fi
    } || {
        echo fail > .status
    }
    """
}

process report_maxbin2_{{ pid }}{

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(tsv) from  intoReport_{{ pid }}

    output:
    {% with task_name="report_maxbin2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_tsv.py"

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