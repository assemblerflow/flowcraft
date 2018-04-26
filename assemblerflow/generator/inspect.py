

class NextflowInspector:

    def __init__(self, trace_file):

        self.trace_file = trace_file
        """
        str: Path to nextflow trace file.
        """

        print(self.trace_file)
