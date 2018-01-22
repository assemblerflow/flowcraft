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

start_channel = """
IN_fastq_raw = Channel.fromFilePairs(params.fastq)
// Channel for expected genome size
IN_genome_size = Channel
                .value(params.genomeSize)
// Channel for minimum coverage threshold
IN_min_coverage = Channel
                .value(params.minCoverage)
                

// Channel for patho_typing
IN_pathoSpecies = Channel
                .value(params.pathoSpecies)

// FASTQC CHANNELS //
// Channel for adapters file
IN_adapters = Channel
                .value(params.adapters)

// TRIMMOMATIC CHANNELS //
IN_trimmomatic_opts = Channel
                .value([params.trimSlidingWindow,
                        params.trimLeading,
                        params.trimTrailing,
                        params.trimMinLength])

// SPADES CHANNELS //
IN_spades_opts = Channel
                .value([params.spadesMinCoverage,
                        params.spadesMinKmerCoverage])
IN_spades_kmers = Channel
                .value(params.spadesKmers)

IN_process_spades_opts = Channel
                .value([params.spadesMinContigLen,
                        params.spadesMinKmerCoverage,
                        params.spadesMaxContigs])

// ASSEMBLY MAPPING CHANNELS //
IN_assembly_mapping_opts = Channel
                .value([params.minAssemblyCoverage,
                        params.AMaxContigs])
    """