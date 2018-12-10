
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Abricate(Process):
    """Abricate mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: None
        - ``ptype``: post_assembly

    It contains one **secondary channel link end**:

        - ``MAIN_assembly`` (alias: ``MAIN_assembly``): Receives the last
        assembly.
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = None

        self.ignore_type = True

        self.status_channels = ["STATUS_abricate", "STATUS_process_abricate"]

        self.params = {
            "abricateDatabases": {
                "default": '["resfinder", "card", "vfdb", "plasmidfinder", '
                           '"virulencefinder", "bacmet"]',
                "description": "Specify the databases for abricate."
            },
            "abricateDataDir": {
                "default": 'null',
                "description": "Specify the full path location of the database "
                               "folders."
            },
            "abricateMinId": {
                "default": '75',
                "description": "Minimum DNA %identity."
            },
            "abricateMinCov": {
                "default": '0',
                "description": "Minimum DNA %coverage."
            }
        }

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})

        self.directives = {
            "abricate": {
                "container": "flowcraft/abricate",
                "version": "0.8.0-3"
            },
            "process_abricate": {
                "container": "flowcraft/abricate",
                "version": "0.8.0-3"
            }
        }


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

        self.directives = {
            "card_rgi": {
                "container": "flowcraft/card_rgi",
                "version": "4.0.2-0.1",
                "memory": "{10.Gb*task.attempt}"
            }
        }

        self.status_channels = [
            "card_rgi"
        ]


class Prokka(Process):
    """Prokka mapping process template interface

    This process is set with:

        - ``input_type``: assembly
        - ``output_type``: None
        - ``ptype``: post_assembly

    It contains one **secondary channel link end**:

        - ``MAIN_assembly`` (alias: ``MAIN_assembly``): Receives the last
        assembly.
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = None

        self.ignore_type = True

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})

        self.params = {
            "centre": {
                "default": "'UMMI'",
                "description": "sequencing centre ID"
            },
            "kingdom": {
                "default": "'Bacteria'",
                "description": "Annotation mode: Archaea|Bacteria|Mitochondria"
                               "|Viruses (default 'Bacteria')"
            },
            "genus": {
                "default": "false",
                "description": "Genus name (default 'Genus'). This also adds"
                               "the --usegenus flag to prokka"
            },
        }

        self.directives = {
            "prokka": {
                "cpus": 2,
                "container": "ummidock/prokka",
                "version": "1.12"
            }
        }


class Diamond(Process):
    """diamond process for protein database queries

        This process is set with:

            - ``input_type``: fasta
            - ``output_type``: None
            - ``ptype``: post_assembly
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = None

        self.params = {
            "pathToDb": {
                "default": 'null',
                "description": "Provide full path for the diamond database. "
                               "If none is provided then will try to fetch from"
                               " the previous process. Default: None"
            },
            "fastaToDb": {
                "default": 'null',
                "description": "Provide the full path for the fasta to "
                               "construct a diamond database. Default: None"
            },
            "blastType": {
                "default": "'blastx'",
                "description": "Defines the type of blast that diamond will do."
                               "Can wither be blastx or blastp. Default: blastx"
            }
        }

        self.directives = {
            "diamond": {
                "container": "flowcraft/diamond",
                "version": "0.9.22-1",
                "memory": "{ 4.GB * task.attempt }",
                "cpus": 2
            }
        }
