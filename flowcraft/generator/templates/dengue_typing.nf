getRef_{{ pid }} = params.reference{{ param_id}} ? "true" : "false"
checkpointReferenceGenome_{{ pid }} = Channel.value(getRef_{{ pid }})

midChan1 = Channel.create()

{{ input_channel }}.join(_LAST_fastq_{{ pid }}).set{ midChan1 }

midChan1.into{ midChan2; midChan3}

midChan3.subscribe { println it }


process dengue_typing_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/dengue_typing/${sample_id}/"


    input:
    set sample_id, file(assembly), file(fastq_pair) from midChan2
    val reference from checkpointReferenceGenome_{{ pid }}

    output:
    file "seq_typing*"
    set sample_id, file("*.fasta") into {{ output_channel}}
    file("*.fa") optional true into _ref_seqTyping_{{ pid }}
    {% with task_name="dengue_typing" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "dengue_typing.py"

}

{{ forks }}

