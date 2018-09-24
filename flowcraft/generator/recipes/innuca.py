try:
    from generator.recipe import Recipe
except ImportError:
    from flowcraft.generator.recipe import Recipe


class Innuca(Recipe):
    """
    Bacterial genome assembly pipeline based on the SPAdes assembler and using
    pre-assembly quality control and read trimming and post-assembly polishing
    with Pilon
    """

    def __init__(self):
        super().__init__()

        # Recipe name
        self.name = "innuca"

        # Recipe pipeline
        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "fastqc " \
                            "check_coverage " \
                            "true_coverage " \
                            "spades " \
                            "process_spades " \
                            "pilon " \
                            "mlst "

        # Recipe parameters and directives
        self.directives = {
            "integrity_coverage": {
                "directives": {"cpus": "1", "memory": "\"2GB\""},
                "params": {"genomeSize": "1", "minCoverage": "15"}
            }
        }
