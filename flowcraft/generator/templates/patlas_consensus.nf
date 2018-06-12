
/**
* A process that creates a consensus from all the outputted json files
*/
process fullConsensus {

    tag { sample_id }

    publishDir 'results/consensus_{{ pid }}/'

    input:
    set sample_id, file(infile_list) from {{ compile_channels }}

    output:
    file "consensus_*.json"

    script:
    template "pATLAS_consensus_json.py"

}