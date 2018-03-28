
process readsDownload_{{ pid }} {

    {% include "post.txt" ignore missing %}

    input:
    file(accession_id) from {{ input_channel }}.splitText()
    file(aspera_key) from Channel.fromPath(params.asperaKey)

    output:
    file "$accession_id/*fastq.gz" into file_list

    {% with task_name="readsDownload" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    echo $accession_id > accession_file.txt
    getSeqENA.py -l ./accession_file.txt -a $aspera_key -o . --SRAopt
--downloadCramBam
    """

}

{{ forks }}
