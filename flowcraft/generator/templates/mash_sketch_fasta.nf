IN_kmerSize_{{ pid }} = Channel.value(params.kmerSize{{ param_id }})
IN_sketchSize_{{ pid }} = Channel.value(params.sketchSize{{ param_id }})

// runs mash sketch
process mashSketchFasta_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(fasta) from {{ input_channel }}
    val kmer from IN_kmerSize_{{ pid }}
    val sketch from IN_sketchSize_{{ pid }}

    output:
    set sample_id, file(fasta) into  {{ output_channel }}
    set sample_id, file("*.msh") into SIDE_mashSketchOutChannel_{{ pid }}
    {% with task_name="mashSketchFasta", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    mash sketch -i -k ${kmer} -s ${sketch} ${fasta}
    """

}

{{ forks }}
