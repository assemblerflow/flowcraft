class ProcessError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ChannelError(Exception):
    def __init__(self, p1, p2, t1, t2):
        self.p1 = p1
        self.p2 = p2
        self.t1 = t1
        self.t2 = t2

    def __str__(self):
        return "The output of the '{}' process ({}) cannot link with the " \
               "input of the '{}' process ({}). Please check the order of " \
               "the processes".format(self.p1, self.p2, self.t1, self.t2)


class SanityError(Exception):
    """
    Class to raise a custom error for sanity checks
    """
    def __init__(self, value):
        self.value = "inSANITY ERROR: " + value

    # def __str__(self):
    #     return repr(self.value)
