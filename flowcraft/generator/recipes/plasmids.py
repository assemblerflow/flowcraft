try:
    from generator.recipe import Recipe
except ImportError:
    from flowcraft.generator.recipe import Recipe


class Plasmids(Recipe):

    def __init__(self):
        super().__init__()

        self.name = "plasmids"

        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "( spades pilon (mash_dist | abricate) |" \
                            "mash_screen | " \
                            "mapping_patlas)"


class PlasmidsMapping(Recipe):

    def __init__(self):
        super().__init__()

        self.name = "plasmids_mapping"

        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "mapping_patlas"


class PlasmidsAssembly(Recipe):

    def __init__(self):
        super().__init__()

        self.name = "plasmids_assembly"

        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "spades " \
                            "pilon " \
                            "mash_dist "


class PlasmidsMash(Recipe):

    def __init__(self):
        super().__init__()

        self.name = "plasmids_mash"

        self.pipeline_str = "integrity_coverage " \
                            "fastqc_trimmomatic " \
                            "mash_screen"
