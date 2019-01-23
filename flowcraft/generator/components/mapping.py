try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Bowtie(Process):
    """
    bowtie2 process to align short paired-end sequencing reads to long reference sequences

        This process is set with:

            - ``input_type``: fastq
            - ``output_type``: bam
            - ``ptype``: mapping

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "bam"

        self.params = {
            "reference": {
                "default": "null",
                "description": "Specifies the reference genome to be provided "
                               "to bowtie2-build."
            },
            "index": {
                "default": "null",
                "description": "Specifies the reference indexes to be provided "
                               "to bowtie2."
            },
            "clearInput": {
                "default": "false",
                "description":
                    "Permanently removes temporary input files. This option "
                    "is only useful to remove temporary files in large "
                    "workflows and prevents nextflow's resume functionality. "
                    "Use with caution."
            }
        }

        self.directives = {
            "bowtie": {
                "container": "flowcraft/bowtie2_samtools",
                "version": "1.0.0-1",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 4
            },
            "bowtie_build": {
                "container": "flowcraft/bowtie2_samtools",
                "version": "1.0.0-1",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 1
            }
        }

        self.status_channels = [
            "bowtie",
            "report_bowtie"
        ]


class RetrieveMapped(Process):
    """
    Samtools process to  to align short paired-end sequencing reads to
    long reference sequences

        This process is set with:

            - ``input_type``: bam
            - ``output_type``: fastq
            - ``ptype``: mapping

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "bam"
        self.output_type = "fastq"

        self.params = {
            "clearInput": {
                "default": "false",
                "description":
                    "Permanently removes temporary input files. This option "
                    "is only useful to remove temporary files in large "
                    "workflows and prevents nextflow's resume functionality. "
                    "Use with caution."
            }
        }

        self.dependencies = ["bowtie"]

        self.directives = {
            "retrieve_mapped": {
                "container": "flowcraft/bowtie2_samtools",
                "version": "1.0.0-1",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 2
            }
        }

        self.status_channels = [
            "retrieve_mapped"
        ]


class InsertSize(Process):
    """
    Determines the sequencing insert size by reads mapping
    to an assembly file

        This process is set with:

            - ``input_type``: fasta
            - ``output_type``: None
            - ``ptype``: mapping

        It contains one **secondary channel link end**:

        - ``MAIN_fq`` (alias: ``_MAIN_assembly``): Receives the FastQ files
        from the last process with ``fastq`` output type.

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = None

        self.link_end.append({"link": "__fastq", "alias": "_LAST_fastq"})

        self.params = {
            "distribution_plot": {
                "default": "false",
                "description": "Produces a distribution plot of the insert sizes."
            },
            "clearInput": {
                "default": "false",
                "description":
                    "Permanently removes temporary input files. This option "
                    "is only useful to remove temporary files in large "
                    "workflows and prevents nextflow's resume functionality. "
                    "Use with caution."
            }
        }

        self.directives = {
            "assembly_mapping_statistics": {
                "container": "flowcraft/bowtie2_samtools",
                "version": "1.0.0-1",
                "memory": "{1.Gb*task.cpus*task.attempt}",
                "cpus": 1
            },
            "insert_size": {
                "container": "flowcraft/plotly",
                "version": "3.5.0-1"
            }
        }

        self.status_channels = [
            "assembly_mapping_statistics",
            "insert_size"
        ]
