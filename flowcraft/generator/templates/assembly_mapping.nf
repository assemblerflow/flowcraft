if ( !params.minAssemblyCoverage{{ param_id }}.toString().isNumber() ){
    if (params.minAssemblyCoverage{{ param_id }}.toString() != 'auto'){
        exit 1, "'minAssemblyCoverage{{ param_id }}' parameter must be a number or 'auto'. Provided value: ${params.minAssemblyCoverage{{ param_id }}}"
    }
}
if ( !params.AMaxContigs{{ param_id }}.toString().isNumber() ){
    exit 1, "'AMaxContigs{{ param_id }}' parameter must be a number. Provide value: '${params.AMaxContigs{{ param_id }}}'"
}

IN_assembly_mapping_opts_{{ pid }} = Channel.value([params.minAssemblyCoverage{{ param_id }},params.AMaxContigs{{ param_id }}])
IN_genome_size_{{ pid }} = Channel.value(params.genomeSize{{ param_id }})


process assembly_mapping_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})

    output:
    set sample_id, file(assembly), 'coverages.tsv', 'coverage_per_bp.tsv', 'sorted.bam', 'sorted.bam.bai' into MAIN_am_out_{{ pid }}
    set sample_id, file("coverage_per_bp.tsv") optional true into SIDE_BpCoverage_{{ pid }}
    {% with task_name="assembly_mapping" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        echo [DEBUG] BUILDING BOWTIE INDEX FOR ASSEMBLY: $assembly >> .command.log 2>&1
        bowtie2-build --threads ${task.cpus} $assembly genome_index >> .command.log 2>&1
        echo [DEBUG] MAPPING READS FROM $fastq >> .command.log 2>&1
        bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -1 ${fastq[0]} -2 ${fastq[1]} -S mapping.sam >> .command.log 2>&1
        echo [DEBUG] CONVERTING AND SORTING SAM TO BAM >> .command.log 2>&1
        samtools sort -o sorted.bam -O bam -@ ${task.cpus} mapping.sam && rm *.sam  >> .command.log 2>&1
        echo [DEBUG] CREATING BAM INDEX >> .command.log 2>&1
        samtools index sorted.bam >> .command.log 2>&1
        echo [DEBUG] ESTIMATING READ DEPTH >> .command.log 2>&1
        parallel -j ${task.cpus} samtools depth -ar {} sorted.bam \\> {}.tab  ::: \$(grep ">" $assembly | cut -c 2-)
        # Insert 0 coverage count in empty files. See Issue #2
        echo [DEBUG] REMOVING EMPTY FILES  >> .command.log 2>&1
        find . -size 0 -print0 | xargs -0 -I{} sh -c 'echo -e 0"\t"0"\t"0 > "{}"'
        echo [DEBUG] COMPILING COVERAGE REPORT  >> .command.log 2>&1
        parallel -j ${task.cpus} echo -n {.} '"\t"' '&&' cut -f3 {} '|' paste -sd+ '|' bc >> coverages.tsv  ::: *.tab
        cat *.tab > coverage_per_bp.tsv
        rm *.tab
        if [ -f "coverages.tsv" ]
        then
            echo pass > .status
        else
            echo fail > .status
        fi
        echo -n "" > .report.json
        echo -n "" > .versions
    } || {
        echo fail > .status
    }
    """
}


/** PROCESS_ASSEMBLY_MAPPING -  MAIN
Processes the results from the assembly_mapping process and filters the
assembly contigs based on coverage and length thresholds.
*/
process process_assembly_mapping_{{ pid }} {

    // Send POST request to platform
    {% with overwrite="false" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(assembly), file(coverage), file(coverage_bp), file(bam_file), file(bam_index) from MAIN_am_out_{{ pid }}
    val opts from IN_assembly_mapping_opts_{{ pid }}
    val gsize from IN_genome_size_{{ pid }}

    output:
    set sample_id, '*_filt.fasta', 'filtered.bam', 'filtered.bam.bai' into {{ output_channel }}
    {% with task_name="process_am" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_assembly_mapping.py"

}

{{ forks }}

