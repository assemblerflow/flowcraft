try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Raxml(Process):
    """mafft to align sequences

            This process is set with:

                - ``input_type``: align
                - ``output_type``: .tree
                - ``ptype``: tree

            """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "align"
        self.output_type = ".tree"

        self.params = {
            "substitutionModel": {
                "default": "'GTRGAMMA'",
                "description": "Substitution model. Option: GTRCAT, GTRCATI, ASC_GTRCAT, GTRGAMMA, ASC_GTRGAMMA etc "
            },
            "seedNumber": {
                "default": "12345",
                "description": "Specify an integer number (random seed) and turn on rapid bootstrapping"
            },
            "bootstrap": {
                "default": "500",
                "description": "Specify the number of alternative runs on distinct starting trees"
            }
        }

        self.directives = {
            "raxml": {
                "container": "flowcraft/raxml",
                "version": "8.2.11-2",
                "cpus": 4,
                "memory": "{ 4.GB * task.attempt }"
            },
            "report_raxml": {
                "container": "flowcraft/raxml",
                "version": "8.2.11-2"
            }
        }

        self.status_channels = [
            "raxml",
            "report_raxml"
        ]


