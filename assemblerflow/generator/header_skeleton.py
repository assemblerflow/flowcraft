header = """#!/usr/bin/nextflow

import Helper
import CheckParams

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

CheckParams.check(params)

nsamples = file(params.fastq).size()
Help.start_info(version, nsamples, "$workflow.start", "$workflow.profile")
    """