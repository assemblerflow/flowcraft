class Help {

    static def start_info(String ver, int fastq, int fasta, String time,
                          String profile) {

        int nsamples = fastq / 2

        println ""
        println "============================================================"
        println "                A S S E M B L E R F L O W"
        println "============================================================"
        println ""
        println " Input FastQ                 : $fastq"
        println " Input Fasta                 : $fasta"
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
        println "                {{ pipeline_name }}"
        println "============================================================"
        println "Built using assemblerflow v{{ version }}"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run {{ nf_file }}"
        println ""
        {% for param, info in help_dict.items() -%}
        println "       --{{ "%-25s" | format(param,) }} {{ info.description }} {{ info.process }}"
        {% endfor %}
    }

}