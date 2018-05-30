try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Mlst(Process):
    """Mlst mapping process template interface

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
        self.output_type = "fasta"

        self.directives = {"mlst": {
            "container": "ummidock/mlst",
        }}

        self.params = {
            "mlstSpecies": {
                "default": "null",
                "description":
                    "Specify the expected species for MLST checking."
            }
        }


class Chewbbaca(Process):
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
            "chewbbaca": {
                "cpus": 4,
                "container": "mickaelsilva/chewbbaca_py3",
                "version": "latest",
            },
            "chewbbaca_batch": {
                "cpus": 4,
                "container": "mickaelsilva/chewbbaca_py3",
                "version": "latest",
            },
            "chewbbacaExtractMLST": {
                "container": "mickaelsilva/chewbbaca_py3",
                "version": "latest"
            }
        }

        self.params = {
            "chewbbacaQueue": {
                "default": "null",
                "description":
                    "Specifiy a queue/partition for chewbbaca. This option"
                    " is only used for grid schedulers. (default: "
                    "$params.chewbbacaQueue)"
            },
            "chewbbacaTraining": {
                "default": "null",
                "description":
                    "Specify the full path to the prodigal training file "
                    "of the corresponding species. (default: "
                    "$params.chewbbacaTraining)"
            },
            "schemaPath": {
                "default": "null",
                "description":
                    "The path to the chewbbaca schema directory. (default: "
                    "$params.schemaPath)"
            },
            "schemaSelectedLoci": {
                "default": "null",
                "description":
                    "The path to the selection of loci in the schema "
                    "directory to be used. If not specified, all loci in the"
                    " schema will be used. (default: "
                    "$params.schemaSelectedLoci)"
            },
            "schemaCore": {
                "default": "null",
                "description": ""
            },
            "chewbbacaJson": {
                "default": "false",
                "description":
                    "If set to True, chewbbaca's allele call output will be "
                    "set to JSON format. (default: $params.chewbbacaJson)"
            },
            "chewbbacaToPhyloviz": {
                "default": "false",
                "description":
                    "If set to True, the ExtractCgMLST module of chewbbaca"
                    " will be executed after the allele calling (default: "
                    "$params.chewbbacaToPhyloviz)",
            },
            "chewbbacaProfilePercentage": {
                "default": 0.95,
                "description":
                    "Specifies the proportion of samples that must be "
                    "present in a locus to save the profile. (default: "
                    "$params.chewbbacaProfilePercentage)"
            },
            "chewbbacaBatch": {
                "default": "false",
                "description":
                    "Specifies whther a chewbbaca run will be performed on the"
                    " complete input batch (all at the same time) or one by "
                    "one."
            }
        }

        self.secondary_inputs = [
            {
                "params": "schemaPath",
                "channel":
                    "if ( !params.schemaPath ){ exit 1, \"'schemaPath' "
                    "parameter missing\"}\n"
                    "if ( params.chewbbacaTraining){"
                    "if (!file(params.chewbbacaTraining).exists()) {"
                    "exit 1, \"'chewbbacaTraining' file was not found: "
                    "'${params.chewbbacaTraining}'\"}}\n"
                    "if ( params.schemaSelectedLoci){"
                    "if (!file(params.schemaSelectedLoci).exists()) {"
                    "exit 1, \"'schemaSelectedLoci' file was not found: "
                    "'${params.schemaSelectedLoci}'\"}}\n"
                    "if ( params.schemaCore){"
                    "if (!file(params.schemaCore).exists()) {"
                    "exit 1, \"'schemaCore' file was not found: "
                    "'${params.schemaCore}'\"}}\n"
                    "IN_schema = Channel.fromPath(params.schemaPath)"
            }
        ]
