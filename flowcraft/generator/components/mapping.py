try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Bowtie(Process):
    """bowtie2 to align short paired-end sequencing reads to long reference sequences

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
    """Samtools process to  to align short paired-end sequencing reads to
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


class Hisat2(Process):
    """HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes (as well as to a single reference genome)

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
                               "to HISAT2."
            },
            "hisat2_index": {
                "default": "null",
                "description": "Specifies the reference indexes to be provided "
                               "to HISAT2."
            },
            "hisat2_index_name": {
                "default": "null",
                "description": "Specifies the reference indexes folder & basename to be provided "
                               "to HISAT2, eg hisat2_index_folder/basename."
            }
        }

        self.directives = {
            "make_hisat2_index": {
                "container": "makaho/hisat2-zstd",
                "version": "latest",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 1
            },
            "hisat2": {
                "container": "makaho/hisat2-zstd",
                "version": "latest",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 4
            }
        }

        self.status_channels = [
            "make_hisat2_index",
            "hisat2"
        ]