class ProcessError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class SanityError(Exception):
    """
    Class to raise a custom error for sanity checks
    """
    def __init__(self, value):
        self.value = "inSANITY ERROR: {}".format(value)

    # def __str__(self):
    #     return repr(self.value)


class InspectionError(Exception):
    def __init__(self, value):
        self.value = "Inspection ERROR: {}".format(value)


class ReportError(Exception):
    def __init__(self, value):
        self.value = "Reports ERROR: {}".format(value)


class RecipeError(Exception):
    def __init__(self, value):
        self.value = "Recipe ERROR: {}".format(value)

    # def __str__(self):
    #     return repr(self.value)

class LogError(Exception):
    def __init__(self, value):
        self.value = "Log ERROR: {}".format(value)
