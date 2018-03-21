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
        self.value = "inSANITY ERROR: " + value

    # def __str__(self):
    #     return repr(self.value)
