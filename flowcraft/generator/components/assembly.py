
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


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

        self.dependencies = ["integrity_coverage"]

        self.params = {
            "spadesMinCoverage": {
                "default": 2,
                "description":
                    "The minimum number of reads to consider an edge in the"
                    " de Bruijn graph during the assembly (default: "
                    "$params.spadesMinCoverage)"
            },
            "spadesMinKmerCoverage": {
                "default": 2,
                "description":
                    "Minimum contigs K-mer coverage. After assembly only "
                    "keep contigs with reported k-mer coverage equal or "
                    "above this value (default: "
                    "$params.spadesMinKmerCoverage)"
            },
            "spadesKmers": {
                "default": "'auto'",
                "description":
                    "If 'auto' the SPAdes k-mer lengths will be determined "
                    "from the maximum read length of each assembly. If "
                    "'default', SPAdes will use the default k-mer lengths. "
                    "(default: $params.spadesKmers)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "spadesOpts",
                "channel":
                    "if ( !params.spadesMinCoverage.toString().isNumber() )"
                    "{ exit 1, \"'spadesMinCoverage' parameter must "
                    "be a number. Provided value: "
                    "${params.spadesMinCoverage}\"}\n"
                    "if ( !params.spadesMinKmerCoverage.toString().isNumber())"
                    "{ exit 1, \"'spadesMinKmerCoverage' parameter must "
                    "be a number. Provided value: "
                    "${params.spadesMinKmerCoverage}\"}\n"
                    "IN_spades_opts = Channel"
                    ".value([params.spadesMinCoverage,"
                    "params.spadesMinKmerCoverage])"
            },
            {
                "params": "spadesKmers",
                "channel":
                    "if ( params.spadesKmers.toString().split(\" \").size() "
                    "<= 1 )"
                    "{ if (params.spadesKmers.toString() != 'auto'){"
                    "exit 1, \"'spadesKmers' parameter must be a sequence "
                    "of space separated numbers or 'auto'. Provided "
                    "value: ${params.spadesKmers}\"} }\n"
                    "IN_spades_kmers = Channel.value(params.spadesKmers)"
            }
        ]

        self.directives = {"spades": {
            "cpus": 4,
            "memory": "{ 5.GB * task.attempt }",
            "container": "ummidock/spades",
            "version": "3.11.1-1",
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
            "container": "ummidock/skesa",
            "version": "0.2.0-3",
            "scratch": "true"
        }}