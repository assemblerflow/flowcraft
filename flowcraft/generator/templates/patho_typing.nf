if ( !params.species{{ param_id }}){ exit 1, "'species' parameter missing" }
if ( params.species{{ param_id }}.toString().split(" ").size() != 2 ){
    exit 1, "'species' parameter must contain two values (e.g.: 'escherichia coli'). Provided value: ${params.species{{ param_id }}}"
}

IN_pathoSpecies_{{ pid }} = Channel.value(params.species{{ param_id }})

process patho_typing_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    errorStrategy "ignore"
    publishDir "results/pathotyping/${sample_id}/"

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}
    val species from IN_pathoSpecies_{{ pid }}

    output:
    file "patho_typing*"
    {% with task_name="patho_typing" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # Prevents read-only issues
        mkdir rematch_temp
        cp -r /NGStools/ReMatCh rematch_temp
        export PATH="\$(pwd)/rematch_temp/ReMatCh:\$PATH"

        patho_typing.py -f \$(pwd)/${fastq_pair[0]} \$(pwd)/${fastq_pair[1]} -o \$(pwd) -j $task.cpus --trueCoverage --species $species
        json_str="{'tableRow':[{'sample':'${sample_id}','data':[{'header':'pathotyping','value':'\$(cat patho_typing.report.txt)','table':'typing'}]}]}"
        echo \$json_str > .report.json

        rm -r rematch_temp
        echo pass > .status

        if [ -s patho_typing.report.txt ];
        then
            echo pass > .status
        else
            echo fail > .status
        fi
    } || {
        echo fail > .status
    }
    """

}

