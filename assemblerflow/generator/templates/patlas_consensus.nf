
/**
* A process that creates a consensus from all the outputted json files
*/
process fullConsensus {

    tag { "Creating consensus json file for: " + id}

    publishDir 'results/consensus/'

    input:
    set id, file(infile_list) from {{ compile_channels }}

    output:
    file "consensus_${id}.json"

    script:
    template "pATLAS_consensus_json.py"

}