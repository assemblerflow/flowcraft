try:
    from generator.recipe import Recipe
except ImportError:
    from flowcraft.generator.recipe import Recipe


class Plasmids(Recipe):
    """
    Plasmid detection pipeline using mapping, mash_screen and assembly with
    SPAdes, with gene annotations with abricate. Outputs json files that
    can be imported into pATLAS.
    """

    def __init__(self):
        super().__init__()

        self.name = "plasmids"

        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "( spades pilon (mash_dist | abricate) |" \
                            "mash_screen | " \
                            "mapping_patlas)"


class PlasmidsMapping(Recipe):
    """
    Plasmid detection pipeline using mapping with bowtie2. Outputs json
    files that can be imported into pATLAS.
    """

    def __init__(self):
        super().__init__()

        self.name = "plasmids_mapping"

        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "mapping_patlas"


class PlasmidsAssembly(Recipe):
    """
    Plasmid detection pipeline using assembly with SPAdes and mash dist.
    Outputs json files that can be imported into pATLAS.
    """

    def __init__(self):
        super().__init__()

        self.name = "plasmids_assembly"

        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "spades " \
                            "pilon " \
                            "mash_dist "


class PlasmidsMash(Recipe):
    """
    Plasmid detection pipeline using mash screen. Outputs json files that can
    be imported into pATLAS.
    """

    def __init__(self):
        super().__init__()

        self.name = "plasmids_mash"

        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "mash_screen"
