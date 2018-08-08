try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Mafft(Process):
    """mafft to align sequences

            This process is set with:

                - ``input_type``: fasta
                - ``output_type``: align
                - ``ptype``: sequence alignment

            """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fasta"
        self.output_type = "align"

        self.params = {
        }

        self.directives = {
            "mafft": {
                "container": "flowcraft/mafft",
                "version": "7.402-1",
                "cpus": 4,
                "memory": "{ 4.GB * task.attempt }"
            }
        }

        self.status_channels = [
            "mafft"
        ]

