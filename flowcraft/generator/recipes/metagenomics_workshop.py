try:
    from generator.recipe import Recipe
except ImportError:
    from flowcraft.generator.recipe import Recipe


class Metagenomics_workshop(Recipe):
    """
    Recipe for ESCMID capacity building workshop 2018 - Clinical Microbiology
    """

    def __init__(self):

        self.name = "metagenomics_workshop"

        self.pipeline_str = "integrity_coverage " \
                            "remove_host " \
                            "fastqc_trimmomatic " \
                            "(megahit " \
                            "assembly_mapping " \
                            "pilon " \
                            "maxbin2 "\
                            "(metaphlan_fa | " \
                            "abricate | " \
                            "mlst ) " \
                            "| metaphlan_fq | " \
                            "metamlst )"

        # Recipe parameters and directives
        self.directives = {
            "abricate": {
                "params": {
                    "abricateDatabases": "\"resfinder\"",
                    "abricateMinCov": '60'
                }
            }
        }
