
process patho_typing {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    errorStrategy "ignore"

    input:
    set fastq_id, file(fastq_pair) from SIDE_PathoType_raw_{{ pid }}
    val species from IN_pathoSpecies

    output:
    file "patho_typing.report.txt"
    set file(".report.json"), file(".status")

    script:
    """
    {
        # Prevents read-only issues
        mkdir rematch_temp
        cp -r /NGStools/ReMatCh rematch_temp
        export PATH="\$(pwd)/rematch_temp/ReMatCh:\$PATH"

        patho_typing.py -f \$(pwd)/${fastq_pair[0]} \$(pwd)/${fastq_pair[1]} -o \$(pwd) -j $task.cpus --trueCoverage --species $species
        json_str="{'typing':{'pathotyping':'\$(cat patho_typing.report.txt)'}}"
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

