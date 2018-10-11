IN_min_contig_lenght_{{ pid }} = Channel.value(params.min_contig_lenght{{ param_id }})
IN_max_iteration_{{ pid }} = Channel.value(params.max_iteration{{ param_id }})
IN_prob_threshold_{{ pid }} = Channel.value(params.prob_threshold{{ param_id }})

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


    output:
    set sample_id, file(assembly), file ('*_maxbin.*.fasta'), file ('bin_status.txt') into binCh_{{ pid }}
    file '*_maxbin.{abundance,log,summary}'
    set sample_id, file("*_maxbin.summary") into intoReport_{{ pid }}
    {% with task_name="maxbin2" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    run_MaxBin.pl -contig ${assembly} -out ${sample_id}_maxbin -reads ${fastq[0]} -reads2 ${fastq[1]} -thread $task.cpus -min_contig_length ${minContigLenght} -max_iteration ${maxIterations} -prob_threshold ${probThreshold}

    if ls *_maxbin.*.fasta 1> /dev/null 2>&1; then echo "true" > bin_status.txt; else echo "false" > false_maxbin.0.fasta; echo "false" > bin_status.txt; fi
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

{{ output_channel }} = Channel.create()

OUT_binned = Channel.create()
OUT_unbinned = Channel.create()

chanA = Channel.create()
chanB = Channel.create()

binCh_{{ pid }}.choice(chanA, chanB){ it -> it[3].text == "false\n" ? 0 : 1 }

chanA.map{ it -> [it[0], it[1]] }.into(OUT_unbinned)

chanB.map{ it -> [it[2].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it[2]]}
    .transpose()
    .map{it -> [it[1].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'),it[1]]}
    .into(OUT_binned)

OUT_binned.mix(OUT_unbinned).into{ {{ output_channel }} }


{{ forks }}