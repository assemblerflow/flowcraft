IN_kmerSize_{{ pid }} = Channel.value(params.kmerSize{{ param_id }})
IN_sketchSize_{{ pid }} = Channel.value(params.sketchSize{{ param_id }})
//IN_genomeSize_{{ pid }} = Channel.value(params.genomeSize{{ param_id }})
IN_minKmer_{{ pid }} = Channel.value(params.minKmer{{ param_id }})


// checks if genomeSize was provided
genomeSize = (params.genomeSize{{ param_id }} == false) ? "" : "-g ${params.genomeSize{{ param_id }}}"

// runs mash sketch
process mashSketchFastq_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(fastq) from {{ input_channel }}
    val kmer from IN_kmerSize_{{ pid }}
    val sketch from IN_sketchSize_{{ pid }}
    val minKmer from IN_minKmer_{{ pid }}

    output:
    set sample_id, file(fastq) into  {{ output_channel }}
    file("*.msh") into SIDE_mashSketchOutChannel_{{ pid }}
    {% with task_name="mashSketchFastq", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    mash sketch -r -k ${kmer} -s ${sketch} -m ${minKmer} ${genomeSize} ${fastq}
    """

}

{{ forks }}