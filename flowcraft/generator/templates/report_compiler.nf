
/** Reports
Compiles the reports from every process
*/
process report {

    tag { sample_id }

    input:
    set sample_id,
            task_name,
            pid,
            report_json,
            version_json,
            trace from {{ compile_channels }}

    output:
    file "*" optional true into master_report

    """
    prepare_reports.py $report_json $version_json $trace $sample_id $task_name 1 $pid $workflow.scriptId $workflow.runName
    """

}

File forkTree = new File("${workflow.projectDir}/.forkTree.json")
File treeDag = new File("${workflow.projectDir}/.treeDag.json")
File js = new File("${workflow.projectDir}/resources/main.js.zip")


forks_channel = forkTree.exists() ?  Channel.fromPath("${workflow.projectDir}/.forkTree.json") : Channel.value(null)
dag_channel = forkTree.exists() ?  Channel.fromPath("${workflow.projectDir}/.treeDag.json") : Channel.value(null)
js_channel = forkTree.exists() ?  Channel.fromPath("${workflow.projectDir}/resources/main.js.zip") : Channel.value(null)

process compile_reports {

    publishDir "pipeline_report/", mode: "copy"

    if ( params.reportHTTP != null ){
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH;"
        afterScript "metadata_POST.sh $params.projectId $params.pipelineId 0 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId 0 \"$params.platformSpecies\""
    }

   input:
   file report from master_report.collect()
   file forks from forks_channel
   file dag from dag_channel
   file js from js_channel

    output:
    file "pipeline_report.json"
    file "pipeline_report.html"
    file "src/main.js"

    script:
    template "compile_reports.py"
}


