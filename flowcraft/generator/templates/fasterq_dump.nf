// check if option file is provided or not
optionFile = (params.option_file{{ param_id }} == false) ? "" :
    "--option-file ${params.option_file{{ param_id }}}"

// process to run fasterq-dump from sra-tools
process fasterqDump_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { accession_id }
    publishDir "reads/${accession_id}/", pattern: "*.fastq*"
    maxRetries 1

    input:
    val accession_id from {{ input_channel }}.splitText(){ it.trim() }.filter{ it.trim() != "" }

    output:
    set accession_id, file("*.fastq*") optional true into {{ output_channel }}
    {% with task_name="fasterqDump", sample_id="accession_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        echo "Downloading the following accession: ${accession_id}"
        fasterq-dump ${accession_id} -e ${task.cpus} -p ${optionFile}
        if [ ${params.compress_fastq{{ param_id }}} = true ]
        then
            echo "Compressing FastQ files..."
            if [ -f ${accession_id}_1.fastq ]
            then
                pigz -p ${task.cpus} ${accession_id}_1.fastq ${accession_id}_2.fastq
            elif [ -f ${accession_id}_3.fastq ]
            then
                echo "No paired end reads were found to compress."
                pigz -p ${task.cpus} ${accession_id}_3.fastq
            else
                echo "FastQ files weren't compressed. Check if FastQ files were downloaded."
            fi
        else
            echo "FastQ files won't be compressed because compress_fastq options was set to: '${params.compress_fastq{{ param_id }}}.'"
        fi
    } || {
        # If exit code other than 0
        if [ \$? -eq 0 ]
        then
            echo "pass" > .status
        else
            echo "fail" > .status
            echo "Could not download accession $accession_id" > .fail
        fi
    }
    """
}
