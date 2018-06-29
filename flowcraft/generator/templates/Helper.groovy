class Help {

    static def start_info(Map info, String time, String profile) {

        println ""
        println "============================================================"
        println "                {{ pipeline_name }}"
        println "============================================================"
        println "Built using flowcraft v{{ version }}"
        println ""
        if (info.containsKey("fastq")){
        int nsamples = info.fastq / 2
        println " Input FastQ                 : $info.fastq"
        println " Input samples               : $nsamples"
        }
        if (info.containsKey("fasta")){
        println " Input Fasta                 : $info.fasta"
        }
        if (info.containsKey("accessions")){
        println " Input accessions            : $info.accessions"
        }
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

    static def print_help(Map params) {

        println ""
        println "============================================================"
        println "                {{ pipeline_name }}"
        println "============================================================"
        println "Built using flowcraft v{{ version }}"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run {{ nf_file }}"
        println ""
        {% for line in help_list -%}
        println "       {{ line }}"
        {% endfor %}
    }

}