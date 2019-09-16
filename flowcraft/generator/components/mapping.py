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
            "retrieve_mapped",
            "renamePE"
        ]


class Bwa(Process):
    """Bwa to align short paired-end sequencing reads to long reference sequences

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
            "bwaIndex": {
                "default": "'s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/BWAIndex/human_g1k_v37_decoy.fasta'",
                "description": "Specifies the reference indexes to be provided "
                               "to bwa."
            }
        }

        self.directives = {
            "bwa": {
                "container": "flowcraft/bwa_samtools",
                "version": "0.7.17-1",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 4
            }
        }

        self.status_channels = [
            "bwa",
        ]


class MarkDuplicates(Process):
    """Identifies duplicate reads.

        This process is set with:

            - ``input_type``: bam
            - ``output_type``: bam
            - ``ptype``: mapping

        """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "bam"
        self.output_type = "bam"

        self.compiler["multiqc"] = ["markDupMultiQC"]

        self.directives = {
            "mark_duplicates": {
                "container": "broadinstitute/gatk",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 4
            }
        }

        self.status_channels = [
            "mark_duplicates"
        ]


class BaseRecalibrator(Process):
    """Detects systematic errors in base quality scores

        This process is set with:

            - ``input_type``: bam
            - ``output_type``: bam
            - ``ptype``: mapping

        """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "bam"
        self.output_type = "bam"

        self.params = {
            "reference": {
                "default": "null",
                "description": "Specifies the name of the FASTA reference genome and index files to be provided "
                               "to BaseRecalibrator."
            },
            "dbsnp": {
                "default": "null",
                "description": "Specifies the dbSNP VCF file to be provided "
                               "to BaseRecalibrator."
            },
            "dbsnpIdx": {
                "default": "null",
                "description": "Specifies the dbSNP VCF index file to be provided "
                               "to BaseRecalibrator."
            },
            "goldenIndel": {
                "default": "null",
                "description": "Specifies the Gold standard INDELs VCF file to be provided "
                               "to BaseRecalibrator."
            },
            "goldenIndelIdx": {
                "default": "null",
                "description": "Specifies the Gold standard INDELs VCF index file to be provided "
                               "to BaseRecalibrator."
            }
        }

        self.directives = {
            "base_recalibrator": {
                "container": "broadinstitute/gatk",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 4
            },
            "apply_bqsr": {
                "container": "broadinstitute/gatk",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 4
            }
        }

        self.status_channels = [
            "base_recalibrator",
            "apply_bqsr"
        ]
