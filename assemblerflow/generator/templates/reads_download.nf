
process readsDownload_{{ pid }} {

    {% include "post.txt" ignore missing %}

    tag { accession_id }

    input:
    val accession_id from {{ input_channel }}.splitText(){ it.trim() }
    each file(aspera_key) from Channel.fromPath(params.asperaKey)

    output:
    set accession_id, file("${accession_id}/*fq.gz") into {{ output_channel }}
    {% with task_name="reads_download", sample_id="accession_id" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    echo "${accession_id}" >> accession_file.txt
    getSeqENA.py -l accession_file.txt -a $aspera_key -o ./ --SRAopt --downloadCramBam
    """

}

{{ forks }}
