class Help {

    static def start_info(String ver, int fastq, String time,
                          String profile) {

        int nsamples = fastq / 2

        println ""
        println "============================================================"
        println " I N N U c a ~ INNUca NextFlow version $ver"
        println "============================================================"
        println ""
        println " Input FastQ                 : $fastq"
        println " Input samples               : $nsamples"
        println " Reports are found in        : ./reports"
        println " Results are found in        : ./results"
        println " Profile                     : $profile"
        println ""
        println "Starting pipeline at $time"
        println ""

    }

    static void complete_info(nextflow.script.WorkflowMetadata wf) {

        println ""
        println "Pipeline execution summary"
        println "=========================="
        println "Completed at                 : $wf.complete"
        println "Duration                     : $wf.duration"
        println "Success                      : $wf.success"
        println "Work directory               : $wf.workDir"
        println "Exit status                  : $wf.exitStatus"
        println ""

    }

    static def print_help(String ver, Map params) {

        println ""
        println "============================================================"
        println " I N N U c a ~ INNUca NextFlow version $ver"
        println "============================================================"
        println "A nextflow implementation of INNUENDO quality control of"
        println "reads, de novo assembly and contigs quality assessment,"
        println "and possible contamination search"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run odiogosilva/innuca-nf"
        println ""
        println "Main options:"
        println "       --fastq                     Patterns for input FastQ files (default: $params.fastq)"
        println "       --genomeSize                Genome size estimate for samples (default: $params.genomeSize)"
        println "       --minCoverage               Minimum coverage (default: $params.minCoverage)"
        println ""
        println "FastQC options:"
        println "       --adapters                  Path to adapters files (default: $params.adapters)"
        println ""
        println "Trimmomatic options:"
        println "       --trimSlidingWindow         Perform sliding window trimming, cutting once the average quality within the window falls below a threshold (default: $params.trimSlidingWindow)"
        println "       --trimLeading               Cut bases off the start of a read, if below a threshold quality (default: $params.trimLeading)"
        println "       --trimTrailing              Cut bases of the end of a read, if below a threshold quality (default: $params.trimTrailing)"
        println "       --trimMinLength             Drop the read if it is below a specified length (default: $params.trimMinLength)"
        println ""
        println "Spades options:"
        println "       --spadesMinCoverage         The minimum number of reads to consider an edge in the de Bruijn graph during the assembly (default: $params.spadesMinCoverage)"
        println "       --spadesMinKmerCoverage     Minimum contigs K-mer coverage. After assembly only keep contigs with reported k-mer coverage equal or above this value (default: $params.spadesMinKmerCoverage)"
        println "       --spadesKmers               If 'auto' the SPAdes k-mer lengths will be determined from the maximum read length of each assembly. If 'default', SPAdes will use the default k-mer lengths. (default: $params.spadesKmers)"
        println "       --spadesMinContigLen        Filter SPAdes contigs for length greater or equal than this value (default: $params.spadesMinContigLen)"
        println "       --spadesMaxContigs          Maximum number of contigs per 1.5 Mb of expected genome size (default: $params.spadesMaxContigs)"
        println ""
        println "Assembly mapping options:"
        println "       --minAssemblyCoverage       In auto, the default minimum coverage for each assembled contig is 1/3 of the assembly mean coverage or 10x, if the mean coverage is below 10x (default: $params.minAssemblyCoverage)"
        println ""
        println "Run-time options:"
        println "       --stopAt                    Specify at which point (process) the pipeline should end (default: params.StopAt)"

    }

}