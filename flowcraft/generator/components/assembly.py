
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Bcalm(Process):
    """Bcalm process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: assembly
        - ``ptype``: assembly

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta"

        self.params = {
            "bcalmKmerSize": {
                "default": 31,
                "description":
                    "size of a kmer"
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

        self.directives = {"bcalm": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "quay.io/biocontainers/bcalm",
            "version": "2.2.0--hd28b015_2",
            "scratch": "true"
        }}


class Spades(Process):
    """Spades process template interface

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
        self.link_start.append("gfa1")

        self.dependencies = ["integrity_coverage"]

        self.params = {
            "spadesMinCoverage": {
                "default": 2,
                "description":
                    "The minimum number of reads to consider an edge in the"
                    " de Bruijn graph during the assembly"
            },
            "spadesMinKmerCoverage": {
                "default": 2,
                "description":
                    "Minimum contigs K-mer coverage. After assembly only "
                    "keep contigs with reported k-mer coverage equal or "
                    "above this value"
            },
            "spadesKmers": {
                "default": "'auto'",
                "description":
                    "If 'auto' the SPAdes k-mer lengths will be determined "
                    "from the maximum read length of each assembly. If "
                    "'default', SPAdes will use the default k-mer lengths. "
            },
            "clearInput": {
                "default": "false",
                "description":
                    "Permanently removes temporary input files. This option "
                    "is only useful to remove temporary files in large "
                    "workflows and prevents nextflow's resume functionality. "
                    "Use with caution."
            },
            "disableRR": {
                "default": "false",
                "description":
                    "disables repeat resolution stage of assembling."
            }
        }

        self.directives = {"spades": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "flowcraft/spades",
            "version": "3.12.0-1",
            "scratch": "true"
        }}


class Skesa(Process):
    """Skesa process template interface
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta"

        self.directives = {"skesa": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "flowcraft/skesa",
            "version": "2.3.0-1",
            "scratch": "true"
        }}

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


class ViralAssembly(Process):
    """
    Process to assemble viral genomes, based on SPAdes and megahit
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta"

        self.dependencies = ["integrity_coverage"]

        self.status_channels = ["va_spades", "va_megahit",
                                "report_viral_assembly"]

        self.link_end.append({"link": "SIDE_max_len", "alias": "SIDE_max_len"})

        self.directives = {"va_spades": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "flowcraft/viral_assembly",
            "version": "0.1-1",
            "scratch": "true"
        }, "va_megahit": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "flowcraft/viral_assembly",
            "version": "0.1-1",
            "scratch": "true"
        }}

        self.params = {
            "minimumContigSize": {
                "default": 10000,
                "description":
                    "Expected genome size in bases"
            },
            "spadesMinCoverage": {
                "default": 2,
                "description":
                    "The minimum number of reads to consider an edge in the"
                    " de Bruijn graph during the assembly"
            },
            "spadesMinKmerCoverage": {
                "default": 2,
                "description":
                    "Minimum contigs K-mer coverage. After assembly only "
                    "keep contigs with reported k-mer coverage equal or "
                    "above this value"
            },
            "spadesKmers": {
                "default": "'auto'",
                "description":
                    "If 'auto' the SPAdes k-mer lengths will be determined "
                    "from the maximum read length of each assembly. If "
                    "'default', SPAdes will use the default k-mer lengths. "
            },
            "megahitKmers": {
                "default": "'auto'",
                "description":
                    "If 'auto' the megahit k-mer lengths will be determined "
                    "from the maximum read length of each assembly. If "
                    "'default', megahit will use the default k-mer lengths. "
                    "(default: $params.megahitKmers)"
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


class Abyss(Process):
    """ABySS process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: assembly
        - ``ptype``: assembly

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta"
        self.link_start.append("gfa1")

        self.params = {
            "abyssKmer": {
                "default": "96",
                "description":
                    "kmer size for assembly."
            }
        }

        self.directives = {"abyss": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "flowcraft/abyss",
            "version": "2.1.1",
            "scratch": "true"
        }}


class Unicycler(Process):
    """Unicycler process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: assembly
        - ``ptype``: assembly

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta"
        self.link_start.append("gfa1")

        self.directives = {"unicycler": {
            "cpus": 4,
            "container": "quay.io/biocontainers/unicycler",
            "version": "0.4.7--py36hdbcaa40_0",
            "scratch": "true"
        }}
