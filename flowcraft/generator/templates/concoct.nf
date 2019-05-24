IN_max_clusters_{{ pid }} = Channel.value(params.clusters{{ param_id }})
IN_length_threshold_{{ pid }} = Channel.value(params.lengthThreshold{{ param_id }})
IN_read_length_{{ pid }} = Channel.value(params.readLength{{ param_id }})
IN_iterations_{{ pid }} = Channel.value(params.iterations{{ param_id }})

clear = params.clearInput{{ param_id }} ? "true" : "false"
checkpointClear_{{ pid }} = Channel.value(clear)

process concoct_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    publishDir "results/assembly/binning/concoct_{{ pid }}/${sample_id}/"

    input:
    set sample_id, file(assembly), file(fastq) from {{ input_channel }}.join(_LAST_fastq_{{ pid }})
    val maxClusters from IN_max_clusters_{{ pid }}
    val read_length from IN_read_length_{{ pid }}
    val length_threshold from IN_length_threshold_{{ pid }}
    val iterations from IN_iterations_{{ pid }}
    val clear from checkpointClear_{{ pid }}

    output:
    set sample_id, file(assembly), file('concoct_output/*.fa') into binCh_{{ pid }}
    set sample_id, file("concoct_output/clustering_merged.csv"), file(assembly) into intoReport_{{ pid }}
    file("concoct_output/*.csv")
    file("concoct_output/*.txt")
    {% with task_name="concoct" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    """
    {
        # cut up the contigs into chunks of 10Kb to mitigate assembly errors and give more weight to larger contigs
        cut_up_fasta.py -c 10000 -o 0 -b ${sample_id}_bedfile -m ${assembly} > ${sample_id}_split_contigs.fasta

        # map reads to cut up assembly
        echo [DEBUG] BUILDING BOWTIE INDEX FOR ASSEMBLY: $assembly >> .command.log 2>&1
        bowtie2-build ${sample_id}_split_contigs.fasta ${sample_id}_split_contigs_index >> .command.log 2>&1
        echo [DEBUG] MAPPING READS FROM $fastq >> .command.log 2>&1
        bowtie2 --threads ${task.cpus} -x ${sample_id}_split_contigs_index -1 ${fastq[0]} -2 ${fastq[1]} -S mapping.sam >> .command.log 2>&1
        echo [DEBUG] CONVERTING AND SORTING SAM TO BAM >> .command.log 2>&1
        samtools sort -o sorted.bam -O bam -@ ${task.cpus} mapping.sam && rm *.sam  >> .command.log 2>&1
        echo [DEBUG] CREATING BAM INDEX >> .command.log 2>&1
        samtools index sorted.bam >> .command.log 2>&1

        # create coverage table for concoct
        concoct_coverage_table.py ${sample_id}_bedfile sorted.bam > ${sample_id}_coverage_file.tab

        # run CONCOCT
        concoct --coverage_file ${sample_id}_coverage_file.tab --composition_file ${sample_id}_split_contigs.fasta \
        -b concoct_output/ -c ${maxClusters} -l ${length_threshold} -r ${read_length } -i ${iterations} -t ${task.cpus}

        # Merge subcontig clustering into original contig clustering
        merge_cutup_clustering.py concoct_output/clustering_*.csv > concoct_output/clustering_merged.csv

        # Extract bins as individual FASTA
        extract_fasta_bins.py --output_path concoct_output/ ${assembly} concoct_output/clustering_merged.csv

        echo pass > .status

        if [ "$clear" = "true" ];
        then
            work_regex=".*/work/.{2}/.{30}/.*"
            file_source1=\$(readlink -f \$(pwd)/${fastq[0]})
            file_source2=\$(readlink -f \$(pwd)/${fastq[1]})
            assembly_file=\$(readlink -f \$(pwd)/${assembly})
            if [[ "\$file_source1" =~ \$work_regex ]]; then
                rm \$file_source1 \$file_source2 \$assembly_file
            fi
        fi
    } || {
        echo fail > .status
    }
    """
}

process report_concoct_{{ pid }}{

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(cluster), file(contigs) from intoReport_{{ pid }}

    output:
    {% with task_name="report_concoct" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "process_concoct.py"

}

// emits one bin per channel
{{ output_channel }} = Channel.create()
binCh_{{ pid }}.map{ it -> [it[2].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it[2]]}
    .transpose()
    .map{it -> [it[1].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'),it[1]]}
    .into({{ output_channel }})

{{ forks }}