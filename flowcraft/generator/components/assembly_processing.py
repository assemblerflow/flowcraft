
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
            "genomeSize": {
                "default": 1,
                "description":
                    "Genome size estimate for the samples in Mb. It is used "
                    "to assess whether an assembly is much larger or smaller "
                    "than expected",
            },
            "skesaMinKmerCoverage": {
                "default": 2,
                "description":
                    "Minimum contigs K-mer coverage. After assembly only keep"
                    " contigs with reported k-mer coverage equal or above "
                    "this value"
            },
            "skesaMinContigLen": {
                "default": 200,
                "description":
                    "Filter contigs for length greater or equal than this "
                    "value"
            },
            "skesaMaxContigs": {
                "default": 100,
                "description":
                    "Maximum number of contigs per 1.5 Mb of expected "
                    "genome size"
            }
        }

        self.directives = {"skesa": {
            "cpus": 1,
            "memory": "'2GB'",
            "container": "flowcraft/skesa",
            "version": "2.1-1",
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
            "genomeSize": {
                "default": 1,
                "description":
                    "Genome size estimate for the samples in Mb. It is used "
                    "to assess whether an assembly is much larger or smaller "
                    "than expected",

            },
            "spadesMinKmerCoverage": {
                "default": 2,
                "description":
                    "Minimum contigs K-mer coverage. After assembly only keep"
                    " contigs with reported k-mer coverage equal or above "
                    "this value"
            },
            "spadesMinContigLen": {
                "default": 200,
                "description":
                    "Filter contigs for length greater or equal than this "
                    "value"
            },
            "spadesMaxContigs": {
                "default": 100,
                "description":
                    "Maximum number of contigs per 1.5 Mb of expected "
                    "genome size"
            }
        }

        self.directives = {"process_spades": {
            "container": "flowcraft/spades",
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
                    " 10x, if the mean coverage is below 10x"
            },
            "AMaxContigs": {
                "default": 100,
                "description":
                    "A warning is issued if the number of contigs is over"
                    "this threshold."
            },
            "genomeSize": {
                "default": 2.1,
                "description":
                    "Genome size estimate for the samples. It is used to "
                    "check the ratio of contig number per genome MB"
            }
        }

        self.directives = {
            "assembly_mapping": {
                "cpus": 4,
                "memory": "{ 5.GB * task.attempt }",
                "container": "flowcraft/bowtie2_samtools",
                "version": "1.0.0-1"
            },
            "process_assembly_mapping": {
                "cpus": 1,
                "memory": "{ 5.GB * task.attempt }",
                "container": "flowcraft/bowtie2_samtools",
                "version": "1.0.0-1"
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
                "container": "flowcraft/pilon",
                "version": "1.22.0-1"
            },
            "pilon_report": {
                "cpus": 1,
                "memory": "{ 7.GB * task.attempt }",
                "container": "flowcraft/pilon",
                "version": "1.22.0-1"
            }
        }
