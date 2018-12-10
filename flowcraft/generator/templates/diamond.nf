// check if any of the parameters it defined before executing the process.
if (!params.pathToDb{{ param_id }} && !params.fastaToDb{{ param_id }})
    exit 1, "'You must specify either a pathToDb or fastaToDb parameter.'"
// checks if both are defined and if so raises an error.
else if (params.pathToDb{{ param_id }} && params.fastaToDb{{ param_id }})
    exit 1, "'Both pathToDb and fastaToDb were given, choose just one.'"

// list of blasts allowed for diamond
allowedBlasts = ["blastp", "blastx"]
// checks if blast type os defined
if (!allowedBlasts.contains(params.blastType{{ param_id }}))
    exit 1, "Provide a valid blast type: blastx or blastp"

process diamond_{{ pid }}  {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir "results/annotation/diamond_{{ pid }}/${sample_id}"

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    file pathToDb from params.pathToDb{{ param_id }} ?
        Channel.fromPath(params.pathToDb{{ param_id }}) : Channel.value("NA")
    file fastaToDb from params.fastaToDb{{ param_id }} ?
        Channel.fromPath(params.fastaToDb{{ param_id }}) : Channel.value("NA")
    val blast from params.blastType{{ param_id }}

    output:
    file "*.txt" into diamondOutputs
    output:
    {% with task_name="diamond"%}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    // Use database when available or otherwise use Fasta file
    if (params.pathToDb{{ param_id }})
        """
        diamond ${blast} -d ${pathToDb} -q ${assembly} \
        -o ${pathToDb}.txt -e 1E-20 -p ${task.cpus} \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend slen sstart send evalue bitscore
        """
    else if (params.fastaToDb{{ param_id }})
        """
        diamond makedb --in ${fastaToDb} -d ${fastaToDb}
        diamond ${blast} -d ${fastaToDb}.dmnd -q ${assembly} \
        -o ${fastaToDb}.txt -e 1E-20 -p ${task.cpus} \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend slen sstart send evalue bitscore
        """

}