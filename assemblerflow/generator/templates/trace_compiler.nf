

process compile_traces {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    input:
    set fastq_id, vals from {{ input_channel }}

   script:
   template "pipeline_status.py"

}

