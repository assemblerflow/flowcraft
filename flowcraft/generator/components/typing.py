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
            "version": "0.1.0-1"
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
        self.output_type = "fasta"

        self.link_start.extend(["_ref_seqTyping"])

        self.params = {
            "reference": {
                "default": "true",
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

        self.status_channels = [
            "dengue_typing"
        ]


class Seqsero2Reads(Process):
    """SeqSero2 for reads process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: None
        - ``ptype``: typing
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = None

        self.directives = {
            "seqsero2_reads": {
                "cpus": 1,
                "memory": "{ 1.GB * task.attempt }",
                "container": "ummidock/seqsero2",
                "version": "alpha-test-1",
                "cache": "false",
                "scratch": "true"
            }
        }


class Seqsero2Assembly(Process):
    """SeqSero2 for assembly process template interface

    This process is set with:

        - ``input_type``: fasta
        - ``output_type``: None
        - ``ptype``: typing
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = None

        self.directives = {
            "seqsero2_assembly": {
                "cpus": 1,
                "memory": "{ 1.GB * task.attempt }",
                "container": "ummidock/seqsero2",
                "version": "alpha-test-1",
                "cache": "false",
                "scratch": "true"
            }
        }


class StxSeqtyping(Process):
    """ecoli_stx_subtyping.py for reads process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: None
        - ``ptype``: typing
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = None

        self.params = {
            "stx2covered": {
                "default": '100',
                "description": "Minimal percentage of sequence covered to consider "
                               "extra stx2 subtypes (value between [0, 100])."
            },
            "stx2identity": {
                "default": '99.5',
                "description": "Minimal sequence identity to consider extra stx2 "
                               "subtypes (value between [0, 100])."
            }
        }

        self.directives = {
            "stx_seqtyping": {
                "cpus": 2,
                "memory": "{ 1.GB * task.cpus * task.attempt }",
                "container": "ummidock/seq_typing",
                "version": "2.2-01",
                "scratch": "true"
            }
        }
