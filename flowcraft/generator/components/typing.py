try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class SeqTyping(Process):
    """

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = None

        self.status_channels = []

        self.link_start = None

        self.directives = {"seq_typing": {
            "cpus": 4,
            "memory": "'4GB'",
            "container": "ummidock/seq_typing",
            "version": "0.1.0-1"
        }}

        self.params = {
            "referenceFileO": {
                "default": "null",
                "description":
                    "Fasta file containing reference sequences. If more"
                    "than one file is passed via the 'referenceFileH parameter"
                    ", a reference sequence for each file will be determined. "
                    "(default: $params.referenceFileO)"
            },
            "referenceFileH": {
                "default": "null",
                "description":
                    "Fasta file containing reference sequences. If more"
                    "than one file is passed via the 'referenceFileO parameter"
                    ", a reference sequence for each file will be determined. "
                    "(default: $params.referenceFileH)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "referenceFileO",
                "channel":
                    "file(params.referenceFileO) ? params.referenceFileO : "
                    "exit(1, \"'referenceFileO' parameter missing\")\n"
                    "IN_refO = Channel"
                    ".fromPath(params.referenceFileO)"
                    "map{ it -> it.exists() ? it : exit(1, \"referenceFileO"
                    " file was not found: '${params.referenceFileO}'\")}"
            },
            {
                "params": "referenceFileH",
                "channel":
                    "file(params.referenceFileH) ? params.referenceFileH : "
                    "exit(1, \"'referenceFileH' parameter missing\")\n"
                    "IN_refH = Channel"
                    ".fromPath(params.referenceFileH)"
                    "map{ it -> it.exists() ? it : exit(1, \"referenceFileH"
                    " file was not found: '${params.referenceFileH}'\")}"
            }
        ]


class PathoTyping(Process):
    """

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = None

        self.ignore_type = True

        self.status_channels = []

        self.params = {
            "species": {
                "default": "null",
                "description":
                    "Species name. Must be the complete species name with"
                    "genus and species, e.g.: 'Yersinia enterocolitica'. "
                    "(default: $params.species)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "species",
                "channel":
                    "if ( !params.species){ exit 1, \"'species' parameter "
                    "missing\" }\n"
                    "if ( params.species.toString().split(\" \").size() != 2 )"
                    "{ exit 1, \"'species' parameter must contain two "
                    "values (e.g.: 'escherichia coli'). Provided value: "
                    "${params.species}\"}\n"
                    "IN_pathoSpecies = Channel.value(params.species)"
            }
        ]

        self.link_start = None
        self.link_end.append({"link": "MAIN_raw",
                              "alias": "SIDE_PathoType_raw"})

        self.directives = {"patho_typing": {
            "cpus": 4,
            "memory": "'4GB'",
            "container": "ummidock/patho_typing",
            "version": "0.3.0-1"
        }}


class Sistr(Process):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = None

        self.directives = {"sistr": {
            "cpus": 4,
            "memory": "'4GB'",
            "container": "ummidock/sistr_cmd",
            "version": "1.0.2"
        }}
