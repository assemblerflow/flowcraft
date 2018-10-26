// Check parameter
if ( !params.bcalmKmerSize{{ param_id }}.toString().isNumber() ){
    exit 1, "'bcalmKmerSize{{ param_id }}' parameter must be a number. Provided value: '${params.bcalmKmes%rSize{{ param_id }}}'"
}

// Clear
clear = params.clearInput{{ param_id }} ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process bcalm_{{ pid }} {
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "reports/assembly/quast_{{pid}}/$sample_id"

    input:
    set sample_id, file(fastq) from {{input_channel}}
    val KmerSize from Channel.value(params.bcalmKmerSize{{param_id}})
    
output:
    file "*.unitig.fa"
    {% with task_name="bcalm" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
	bcalm -in $fastq -out unitig -kmer-size $KmerSize"

  	if [ "$clear" = "true" ];
	then
    	    find . -type f  -print | egrep "work/.*(h5)|(glue)" | xargs -L 1 rm
	fi
    }
    """
}

