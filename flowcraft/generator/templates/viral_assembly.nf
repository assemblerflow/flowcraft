//MAIN INPUT - FASTQ FILES
spades_in = Channel.create()
megahit_in = Channel.create()
{{ input_channel }}.into{ spades_in; megahit_in }

//EXPECTED GENOME SIZE
if ( !params.minimumContigSize{{ param_id }}.toString().isNumber() ){
    exit 1, "'minimumContigSize{{ param_id }}' parameter must be a number. Provided value: '${params.minimumContigSize{{ param_id }}}'"
}

//SPADES OPTIONS
if ( !params.spadesMinCoverage{{ param_id }}.toString().isNumber() ){
    exit 1, "'spadesMinCoverage{{ param_id }}' parameter must be a number. Provided value: '${params.spadesMinCoverage{{ param_id }}}'"
}
if ( !params.spadesMinKmerCoverage{{ param_id }}.toString().isNumber()){
    exit 1, "'spadesMinKmerCoverage{{ param_id }}' parameter must be a number. Provided value: '${params.spadesMinKmerCoverage{{ param_id }}}'"
}

if ( params.spadesKmers{{ param_id }}.toString().split(" ").size() <= 1 ){
    if (params.spadesKmers{{ param_id }}.toString() != 'auto'){
        exit 1, "'spadesKmers{{ param_id }}' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.spadesKmers{{ param_id }}}"
    }
}

clear = params.clearAtCheckpoint ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

//MEGAHIT OPTIONS
if ( params.megahitKmers{{ param_id }}.toString().split(" ").size() <= 1 ){
    if (params.megahitKmers{{ param_id }}.toString() != 'auto'){
        exit 1, "'megahitKmers{{ param_id }}' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.megahitKmers{{ param_id }}}"
    }
}

//SPADES INPUT CHANNELS
IN_spades_opts_{{ pid }} = Channel.value([params.spadesMinCoverage{{ param_id }},params.spadesMinKmerCoverage{{ param_id }}])
IN_spades_kmers_{{ pid }} = Channel.value(params.spadesKmers{{ param_id }})

//MEGAGIT INPUT CHANNELS
IN_megahit_kmers_{{ pid }} = Channel.value(params.megahitKmers{{ param_id }})

SIDE_max_len_spades = Channel.create()
SIDE_max_len_megahit = Channel.create()
SIDE_max_len_{{ pid }}.into{SIDE_max_len_spades ; SIDE_max_len_megahit}

process va_spades_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    validExitStatus 0,1

    tag { sample_id }
    publishDir 'results/assembly/spades_{{ pid }}/', pattern: '*_spades*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from spades_in.join(SIDE_max_len_spades)
    val opts from IN_spades_opts_{{ pid }}
    val kmers from IN_spades_kmers_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, file({task.exitStatus == 1 ? ".exitcode" : '*_spades*.fasta'}) into assembly_spades
    {% with task_name="va_spades" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "spades.py"

}

class VerifyCompletness {

    public static boolean contigs(String filename, int threshold){
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        boolean result = processContigs(reader, threshold);
        reader.close()

        return result;
    }

    private static boolean processContigs(BufferedReader reader, int threshold){
        String line;
        int lineThreshold = 0;
        List splittedLine

        while ((line = reader.readLine()) != null) {
            if (line.startsWith('>')) {
                splittedLine = line.split('_')
                lineThreshold = splittedLine[3].toInteger()
                if(lineThreshold >= threshold) {
                    return true;
                }
             }
        }

        return false;
    }
}

megahit = Channel.create()
good_assembly = Channel.create()
assembly_spades.choice(good_assembly, megahit){a -> a[1].toString() == "null" ? false : VerifyCompletness.contigs(a[1].toString(), params.minimumContigSize{{ param_id }}.toInteger()) == true ? 0 : 1}


process va_megahit_{{ pid }}  {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    publishDir 'results/assembly/megahit_{{ pid }}/', pattern: '*_megahit*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from megahit_in.join(megahit).map{ ot -> [ot[0], ot[1]] }.join(SIDE_max_len_megahit)
    val kmers from IN_megahit_kmers_{{ pid }}

    output:
    set sample_id, file('*megahit*.fasta') into megahit_assembly
    {% with task_name="va_megahit" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "megahit.py"

}


good_assembly.mix(megahit_assembly).set{ {{ output_channel }} }


{{ forks }}