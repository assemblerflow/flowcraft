IN_fragLen_{{ pid }} = Channel.value(params.fragLen{{ param_id }})

// runs fast ani for multiple comparisons (many to many mode)
process fastAniMatrix_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { sample_id }

     publishDir 'results/fast_ani/fast_ani_{{ pid }}/',

    input:
    set sample_id, file(fasta) from {{ input_channel }}
    val fragLenValue from IN_fragLen_{{ pid }}

    output:
    set sample_id, fasta, file("*.out")
    {% with task_name="fastAniMatrix", sample_id="sample_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    mkdir fasta_store
    fasta_spliter.py ${fasta}
    fastANI --ql files_fastani.txt --rl files_fastani.txt \
    -t ${task.cpus} --fragLen ${fragLenValue} \
    -o ${sample_id.take(sample_id.lastIndexOf("."))}_fastani.out
    """

}
