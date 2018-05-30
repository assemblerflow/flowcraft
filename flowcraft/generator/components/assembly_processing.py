
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class ProcessSkesa(Process):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = "fasta"

        self.params = {
            "skesaMinKmerCoverage": {
                "default": 2,
                "description":
                    "Minimum contigs K-mer coverage. After assembly only keep"
                    " contigs with reported k-mer coverage equal or above "
                    "this value (default: $params.skesaMinKmerCoverage)"
            },
            "skesaMinContigLen": {
                "default": 200,
                "description":
                    "Filter contigs for length greater or equal than this "
                    "value (default: $params.skesaMinContigLen)"
            },
            "skesaMaxContigs": {
                "default": 100,
                "description":
                    "Maximum number of contigs per 1.5 Mb of expected "
                    "genome size (default: $params.skesaMaxContigs)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "processSkesaOpts",
                "channel":
                    "if ( !params.skesaMinKmerCoverage.toString().isNumber() )"
                    "{ exit 1, \"'skesaMinKmerCoverage' parameter must "
                    "be a number. Provided value: "
                    "${params.skesaMinKmerCoverage}\"}\n"
                    "if ( !params.skesaMinContigLen.toString().isNumber() )"
                    "{ exit 1, \"'skesaMinContigLen' parameter must "
                    "be a number. Provided value: "
                    "${params.skesaMinContigLen}\"}\n"
                    "if ( !params.skesaMaxContigs.toString().isNumber() )"
                    "{ exit 1, \"'skesaMaxContigs' parameter must "
                    "be a number. Provided value: "
                    "${params.skesaMaxContigs}\"}\n"
                    "IN_process_skesa_opts = Channel"
                    ".value([params.skesaMinContigLen,"
                    "params.skesaMinKmerCoverage,params.skesaMaxContigs])"
            }
        ]

        self.directives = {"skesa": {
            "cpus": 1,
            "memory": "'2GB'",
            "container": "ummidock/skesa",
            "version": "0.2.0-3",
        }}


class ProcessSpades(Process):
    """Process spades process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: assembly
        - ``ptype``: post_assembly

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = "fasta"

        self.params = {
            "spadesMinKmerCoverage": {
                "default": 2,
                "description":
                    "Minimum contigs K-mer coverage. After assembly only keep"
                    " contigs with reported k-mer coverage equal or above "
                    "this value (default: $params.spadesMinKmerCoverage)"
            },
            "spadesMinContigLen": {
                "default": 200,
                "description":
                    "Filter contigs for length greater or equal than this "
                    "value (default: $params.spadesMinContigLen)"
            },
            "spadesMaxContigs": {
                "default": 100,
                "description":
                    "Maximum number of contigs per 1.5 Mb of expected "
                    "genome size (default: $params.spadesMaxContigs)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "processSpadesOpts",
                "channel":
                    "if ( !params.spadesMinKmerCoverage.toString().isNumber())"
                    "{ exit 1, \"'spadesMinKmerCoverage' parameter must "
                    "be a number. Provided value: "
                    "${params.spadesMinKmerCoverage}\"}\n"
                    "if ( !params.spadesMinContigLen.toString().isNumber() )"
                    "{ exit 1, \"'spadesMinContigLen' parameter must "
                    "be a number. Provided value: "
                    "${params.spadesMinContigLen}\"}\n"
                    "if ( !params.spadesMaxContigs.toString().isNumber() )"
                    "{ exit 1, \"'spadesMaxContigs' parameter must "
                    "be a number. Provided value: "
                    "${params.spadesMaxContigs}\"}\n"
                    "IN_process_spades_opts = Channel"
                    ".value([params.spadesMinContigLen, "
                    "params.spadesMinKmerCoverage, params.spadesMaxContigs])"
            }
        ]

        self.directives = {"process_spades": {
            "container": "ummidock/spades",
            "version": "3.11.1-1"
        }}


class AssemblyMapping(Process):
    """Assembly mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: assembly
        - ``ptype``: post_assembly

    It contains one **secondary channel link end**:

        - ``MAIN_fq`` (alias: ``_MAIN_assembly``): Receives the FastQ files
        from the last process with ``fastq`` output type.

    It contains two **status channels**:

        - ``STATUS_am``: Status for the assembly_mapping process
        - ``STATUS_amp``: Status for the process_assembly_mapping process
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = "fasta"

        self.status_channels = ["STATUS_assembly_mapping",
                                "STATUS_process_am"]

        self.link_start.append("SIDE_BpCoverage")
        self.link_end.append({"link": "__fastq", "alias": "_LAST_fastq"})

        self.params = {
            "minAssemblyCoverage": {
                "default": "'auto'",
                "description":
                    "In auto, the default minimum coverage for each "
                    "assembled contig is 1/3 of the assembly mean coverage or"
                    " 10x, if the mean coverage is below 10x (default: "
                    "$params.minAssemblyCoverage)"
            },
            "AMaxContigs": {
                "default": 100,
                "description":
                    "A warning is issues if the number of contigs is over"
                    "this threshold"
            },
            "genomeSize": {
                "default": 2.1,
                "description":
                    "Genome size estimate for the samples. It is used to "
                    "check the ratio of contig number per genome MB "
                    "(default: $params.genomeSize)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "assemblyMappingOpts",
                "channel":
                    "if ( !params.minAssemblyCoverage.toString().isNumber() )"
                    "{ if (params.minAssemblyCoverage.toString() != 'auto'){"
                    "exit 1, \"'minAssemblyCoverage' parameter must be a"
                    " number or 'auto'. Provided value: "
                    "${params.minAssemblyCoverage}\"} }\n"
                    "if ( !params.AMaxContigs.toString().isNumber() )"
                    "{ exit 1, \"'AMaxContigs' parameter must be a number."
                    "Provide value: '${params.AMaxContigs}'\"}\n"
                    "IN_assembly_mapping_opts = Channel"
                    ".value([params.minAssemblyCoverage,params.AMaxContigs])"
            },
            {
                "params": "genomeSize",
                "channel":
                    "if ( !params.genomeSize.toString().isNumber() )"
                    "{ exit 1, \"'genomeSize' parameter must be a number."
                    "Provide value: '${params.genomeSize}'\"}\n"
                    "IN_genome_size = Channel.value(params.genomeSize)"
            }
        ]

        self.directives = {
            "assembly_mapping": {
                "cpus": 4,
                "memory": "{ 5.GB * task.attempt }",
                "container": "ummidock/bowtie2_samtools",
                "version": "1.0.0-2"
            },
            "process_assembly_mapping": {
                "cpus": 1,
                "memory": "{ 5.GB * task.attempt }",
                "container": "ummidock/bowtie2_samtools",
                "version": "1.0.0-2"
            }
        }


class Pilon(Process):
    """Pilon mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: assembly
        - ``ptype``: post_assembly

    It contains one **dependency process**:

        - ``assembly_mapping``: Requires the BAM file generated by the
        assembly mapping process
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = "fasta"

        self.dependencies = ["assembly_mapping"]
        self.status_channels = ["STATUS_pilon", "STATUS_pilon_report"]

        self.link_end.append({"link": "SIDE_BpCoverage",
                              "alias": "SIDE_BpCoverage"})

        self.directives = {
            "pilon": {
                "cpus": 4,
                "memory": "{ 7.GB * task.attempt }",
                "container": "ummidock/pilon",
                "version": "1.22.0-2"
            },
            "pilon_report": {
                "cpus": 1,
                "memory": "{ 7.GB * task.attempt }",
                "container": "ummidock/pilon",
                "version": "1.22.0-2"
            }
        }
