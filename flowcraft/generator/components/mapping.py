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

class Retrieve_mapped(Process):
    """Samtools process to  to align short paired-end sequencing reads to long reference sequences

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
