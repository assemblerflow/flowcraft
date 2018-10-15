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

        self.link_start = None

        self.directives = {"seq_typing": {
            "cpus": 4,
            "memory": "'4GB'",
            "container": "flowcraft/seq_typing",
            "version": "2.0-1"
        }}

        self.params = {
            "referenceFileO": {
                "default": "null",
                "description":
                    "Fasta file containing reference sequences. If more"
                    "than one file is passed via the 'referenceFileH parameter"
                    ", a reference sequence for each file will be determined. "
            },
            "referenceFileH": {
                "default": "null",
                "description":
                    "Fasta file containing reference sequences. If more"
                    "than one file is passed via the 'referenceFileO parameter"
                    ", a reference sequence for each file will be determined. "
            }
        }


class PathoTyping(Process):
    """

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = None

        self.ignore_type = True

        self.params = {
            "species": {
                "default": "null",
                "description":
                    "Species name. Must be the complete species name with"
                    "genus and species, e.g.: 'Yersinia enterocolitica'. "
            }
        }

        self.link_start = None
        self.link_end.append({"link": "MAIN_raw",
                              "alias": "SIDE_PathoType_raw"})

        self.directives = {"patho_typing": {
            "cpus": 4,
            "memory": "'4GB'",
            "container": "flowcraft/patho_typing",
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


class Momps(Process):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = None

        self.link_end.append({"link": "__fastq", "alias": "_LAST_fastq"})

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

        self.directives = {
            "momps": {
                "cpus": 3,
                "memory": "'4GB'",
                "container": "flowcraft/momps",
                "version": "0.1.1-1"
            }
        }


class DengueTyping(Process):
    """

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = None

        self.link_start = None

        self.params = {
            "reference": {
                "default": "false",
                "description":
                    "Retrieves the sequence of the closest reference."
            }
        }

        self.directives = {"dengue_typing": {
            "cpus": 4,
            "memory": "'4GB'",
            "container": "flowcraft/seq_typing",
            "version": "2.0-1"
        }}


