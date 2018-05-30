
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

        self.status_channels = ["STATUS_abricate"]

        self.params = {
            "abricateDatabases": {
                "default": '["resfinder", "card", "vfdb", "plasmidfinder", '
                           '"virulencefinder"]',
                "description": "Specify the databases for abricate."
            }
        }

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})

        self.directives = {
            "abricate": {
                "container": "ummidock/abricate",
                "version": "0.8.0-1"
            },
            "process_abricate": {
                "container": "ummidock/abricate",
                "version": "0.8.0-1"
            }
        }


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

        self.directives = {
            "prokka": {
                "cpus": 2,
                "container": "ummidock/prokka-nf",
                "version": "1.12.0-2"
            }
        }
