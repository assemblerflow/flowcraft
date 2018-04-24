
try:
    from generator.process import Process
except ImportError:
    from assemblerflow.generator.process import Process

class remove_host(Process):

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
                "memory": "{10.Gb*task.attempt}",
                "cpus": 3
            }
        }

        self.status_channels = [
            "remove_host"
        ]

class metaspades(Process):
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


class card_rgi(Process):

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
                "memory": "{10.Gb*task.attempt}",
            }
        }

        self.status_channels = [
            "card_rgi"
        ]



