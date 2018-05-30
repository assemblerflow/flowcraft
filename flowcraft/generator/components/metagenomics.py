
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class CardRgi(Process):
    """card's rgi process template interface

        This process is set with:

            - ``input_type``: fasta
            - ``output_type``: txt
            - ``ptype``: resistance gene detection (assembly)
        """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = "txt"

        self.params = {
            "alignmentTool": {
                "default": "'DIAMOND'",
                "description": "Specifies the alignment tool."
                               "Options: DIAMOND or BLAST"
            }
        }

        self.secondary_inputs = [
            {
                "params": "alignmentTool",
                "channel": "IN_alignment_tool = Channel.value(params.alignmentTool)"
            }
        ]

        self.directives = {
            "card_rgi": {
                "container": "cimendes/card_rgi",
                "version": "4.0.2",
                "memory": "{10.Gb*task.attempt}"
            }
        }

        self.status_channels = [
            "card_rgi"
        ]


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
                "default": "'minikraken_20141208'",
                "description": "Specifies kraken database."
            }
        }

        self.secondary_inputs = [
            {
                "params": "krakenDB",
                "channel": "IN_kraken_DB = Channel.value(params.krakenDB)"
            }
        ]

        self.directives = {
            "kraken": {
                "container": "ummidock/kraken_metagenomics",
                "version": "v0.10.6-unreleased",
                "memory": "{2.Gb*task.attempt}",
                "cpus": 3
            }
        }

        self.status_channels = [
            "kraken"
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
            }
        }

        self.secondary_inputs = [
            {
                "params": "megahitKmers",
                "channel":
                    "if ( params.megahitKmers.toString().split(\" \").size() "
                    "<= 1 )"
                    "{ if (params.megahitKmers.toString() != 'auto'){"
                    "exit 1, \"'megahitKmers' parameter must be a sequence "
                    "of space separated numbers or 'auto'. Provided "
                    "value: ${params.megahitKmers}\"} }\n"
                    "IN_megahit_kmers = Channel.value(params.megahitKmers)"
            }
        ]

        self.directives = {"megahit": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "cimendes/megahit",
            "version": "v1.1.3-0.1",
            "scratch": "true"
        }}


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
            }
        }

        self.secondary_inputs = [
            {
                "params": "metaspadesKmers",
                "channel":
                    "if ( params.metaspadesKmers.toString().split(\" \").size() "
                    "<= 1 )"
                    "{ if (params.metaspadesKmers.toString() != 'auto'){"
                    "exit 1, \"'metaspadesKmers' parameter must be a sequence "
                    "of space separated numbers or 'auto'. Provided "
                    "value: ${params.metaspadesKmers}\"} }\n"
                    "IN_metaspades_kmers = Channel.value(params.metaspadesKmers)"
            }
        ]

        self.directives = {"metaspades": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "ummidock/spades",
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
                "default": "'/MidasDB/midas_db_v1.2'",
                "description": "Specifies Midas database."
            }
        }

        self.secondary_inputs = [
            {
                "params": "midasDB",
                "channel": "IN_midas_DB = Channel.value(params.midasDB)"
            }
        ]

        self.directives = {
            "midas_species": {
                "container": "cimendes/midas",
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
            }
        }

        self.secondary_inputs = [
            {
                "params": "refIndex",
                "channel": "IN_index_files = Channel.value(params.refIndex)"
            }
        ]

        self.directives = {
            "remove_host": {
                "container": "cimendes/bowtie2_hg19",
                "version": "2.0",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 3
            }
        }

        self.status_channels = [
            "remove_host"
        ]

class MetaProb(Process):

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
            "mashScreen": {
                "container": "",
                "version": "1.0.0",
                "cpus": 1,
                "memory": "{ 4.GB * task.attempt }"
            },
            "mashOutputJson": {
                "container": "",
                "version": "1.0.0",
                "cpus": 1,
                "memory": "'4GB'"
            }
        }

        self.status_channels = [
            "metaProb"
        ]


