try:
    from generator.recipe import Recipe
except ImportError:
    from flowcraft.generator.recipe import Recipe


class Denim(Recipe):

    def __init__(self):

        self.name = "denim"
