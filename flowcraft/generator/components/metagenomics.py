
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Kraken(Process):
    """kraken process template interface

            This process is set with:

                - ``input_type``: fastq
                - ``output_type``: txt
                - ``ptype``: taxonomic classification
    """
    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "txt"

        self.params = {
            "krakenDB": {
                "default": "'minikraken_20171013_4GB'",
                "description": "Specifies kraken database."
            }
        }

        self.directives = {
            "kraken": {
                "container": "flowcraft/kraken",
                "version": "1.0-0.1",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 3
            }
        }

        self.status_channels = [
            "kraken"
        ]

class Kraken2(Process):
    """kraken2 process template interface

            This process is set with:

                - ``input_type``: fastq
                - ``output_type``: txt
                - ``ptype``: taxonomic classification
    """
    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = None

        self.params = {
            "kraken2DB": {
                "default": "'minikraken2_v1_8GB'",
                "description": "Specifies kraken2 database. Requires full path if database not on "
                               "KRAKEN2_DB_PATH."
            }
        }

        self.directives = {
            "kraken2": {
                "container": "flowcraft/kraken2",
                "version": "2.0.7-1",
                "memory": "{8.Gb*task.attempt}",
                "cpus": 4
            }
        }

        self.status_channels = [
            "kraken2"
        ]


class Maxbin2(Process):
    """MaxBin2, a metagenomics binning software

            This process is set with:

                - ``input_type``: assembly
                - ``output_type``: assembly
                - ``ptype``: post_assembly

            It contains one **secondary channel link end**:

                - ``MAIN_fq`` (alias: ``_MAIN_assembly``): Receives the FastQ files
                from the last process with ``fastq`` output type.

            """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = "fasta"

        self.link_end.append({"link": "__fastq", "alias": "_LAST_fastq"})

        self.params = {
            "min_contig_lenght": {
                "default": 1000,
                "description": "minimum contig length. Default: 1000"
            },
            "max_iteration": {
                "default": 50,
                "description": "maximum Expectation-Maximization algorithm"
                               "iteration number. Default: 50"
            },
            "prob_threshold": {
                "default": 0.9,
                "description": "probability threshold for EM final classification."
                               "Default: 0.9"
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
            "maxbin2": {
                "container": "flowcraft/maxbin2",
                "version": "2.2.4-1",
                "cpus": 3,
                "memory": "{ 5.GB * task.attempt }"
            }
        }

        self.status_channels = [
            "maxbin2",
            "report_maxbin2"
        ]


class Megahit(Process):
    """megahit process template interface

        This process is set with:

            - ``input_type``: fastq
            - ``output_type``: assembly
            - ``ptype``: assembly

        It contains one **secondary channel link end**:

            - ``SIDE_max_len`` (alias: ``SIDE_max_len``): Receives max read length
        """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta"

        self.link_end.append({"link": "SIDE_max_len", "alias": "SIDE_max_len"})

        self.dependencies = ["integrity_coverage"]

        self.params = {
            "megahitKmers": {
                "default": "'auto'",
                "description":
                    "If 'auto' the megahit k-mer lengths will be determined "
                    "from the maximum read length of each assembly. If "
                    "'default', megahit will use the default k-mer lengths. "
                    "(default: $params.megahitKmers)"
            },
            "fastg": {
                "default": "false",
                "description":
                    "Converts megahit intermediate contigs to fastg"

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

        self.directives = {"megahit": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "flowcraft/megahit",
            "version": "1.1.3-0.1",
            "scratch": "true"
        },
            "megahit_fastg": {
                "container": "flowcraft/megahit",
                "version": "1.1.3-0.1",
            }
        }

        self.status_channels = [
            "megahit",
            "megahit_fastg"
        ]


class Metaspades(Process):
    """Metaspades process template interface

        This process is set with:

            - ``input_type``: fastq
            - ``output_type``: assembly
            - ``ptype``: assembly

        It contains one **secondary channel link end**:

            - ``SIDE_max_len`` (alias: ``SIDE_max_len``): Receives max read length
        """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta"

        self.link_end.append({"link": "SIDE_max_len", "alias": "SIDE_max_len"})

        self.dependencies = ["integrity_coverage"]

        self.params = {
            "metaspadesKmers": {
                "default": "'auto'",
                "description":
                    "If 'auto' the metaSPAdes k-mer lengths will be determined "
                    "from the maximum read length of each assembly. If "
                    "'default', metaSPAdes will use the default k-mer lengths. "
                    "(default: $params.metaspadesKmers)"
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

        self.directives = {"metaspades": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "flowcraft/spades",
            "version": "3.11.1-1",
            "scratch": "true"
        }}


class Midas_species(Process):
    """Midas species process template interface

            This process is set with:

                - ``input_type``: fastq
                - ``output_type``: txt
                - ``ptype``: taxonomic classification (species)
    """
    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "txt"

        self.params = {
            "midasDB": {
                "default": "null",
                "description": "Specifies Midas database."
            }
        }

        self.directives = {
            "midas_species": {
                "container": "flowcraft/midas",
                "version": "1.3.2-0.1",
                "memory": "{2.Gb*task.attempt}",
                "cpus": 3
            }
        }

        self.status_channels = [
            "midas_species"
        ]


class RemoveHost(Process):
    """bowtie2 to remove host reads process template interface

        This process is set with:

            - ``input_type``: fastq
            - ``output_type``: fastq
            - ``ptype``: removal os host reads

        """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.params = {
            "refIndex": {
                "default": "'/index_hg19/hg19'",
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
            "remove_host": {
                "container": "flowcraft/remove_host",
                "version": "2-0.1",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 3
            }
        }

        self.status_channels = [
            "remove_host",
            "report_remove_host"
        ]

class Metaprob(Process):
    """MetaProb to bin metagenomic reads interface

            This process is set with:

                - ``input_type``: fastq
                - ``output_type``: csv
                - ``ptype``: binning of reads

            """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "csv"

        self.params = {
            "feature": {
                "default": 1,
                "description": "Feature used to compute. Default: 1"
            },
            "metaProbQMer": {
                "default": 5,
                "description": "Threshold of shared q-mer to create graph "
                               "adiacences. Default: 5"
            }
        }

        self.directives = {
            "metaProb": {
                "container": "flowcraft/metaprob",
                "version": "2-1",
                "cpus": 1,
                "memory": "{ 30.GB * task.attempt }"
            }
        }

        self.status_channels = [
            "metaProb"
        ]


class SplitAssembly(Process):
    """Component to filter metagenomic assemblies by contig size
    If the contig is larger than $param.size, it gets separated
    from the original assembly to continue the processes downstream
    of the pipeline.

            This process is set with:

                - ``input_type``: fasta
                - ``output_type``: fasta
                - ``ptype``: assembly filter

            """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = "fasta"


        self.params = {
            "size": {
                "default": "null",
                "description": "Minimum contig size"
            }
        }


        self.directives = {
            "split_assembly": {
                "cpus": 1,
                "memory": "{ 1.GB * task.attempt }"
            }
        }

        self.status_channels = [
            "split_assembly"
        ]

