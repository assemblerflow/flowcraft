IN_adapter_{{ pid }} = Channel.value(params.adapter{{ param_id }})


process filter_poly_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    errorStrategy { task.exitStatus == 120 ? 'ignore' : 'retry' }

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val adapter from IN_adapter_{{ pid }}

    output:
    set sample_id , file("${sample_id}_filtered_{1,2}.fastq.gz") into {{ output_channel }}
    {% with task_name="filter_poly" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    gunzip -c ${fastq_pair[0]} >  ${sample_id}_1.fq
    gunzip -c ${fastq_pair[1]} >  ${sample_id}_2.fq

    for seqfile in *.fq;
    do if [ ! -s \$seqfile  ]
    then
        echo \$seqfile is empty && exit 120
    fi
    done

    prinseq-lite.pl --fastq ${sample_id}_1.fq  --fastq2 ${sample_id}_2.fq  --custom_params "${adapter}" -out_format 3 -out_good ${sample_id}_filtered

    gzip ${sample_id}_filtered_*.fastq

    """
}

{{ forks }}

