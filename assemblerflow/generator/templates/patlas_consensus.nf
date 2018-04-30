
/**
* A process that creates a consensus from all the outputted json files
*/
process fullConsensus {

    tag { "Creating consensus json file for: " + sample_id}

    publishDir 'results/consensus/'

    input:
    set sample_id, file(infile_list) from {{ compile_channels }}

    output:
    file "consensus_${sample_id}.json"

    script:
    template "pATLAS_consensus_json.py"

}