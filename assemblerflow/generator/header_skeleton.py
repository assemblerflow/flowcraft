header = """#!/usr/bin/env nextflow

import Helper

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(version, params)
    exit 0
}

nsamples = file(params.fastq).size()
nfasta = file(params.fasta).size()
Help.start_info(version, nsamples, nfasta, "$workflow.start", "$workflow.profile")
    """