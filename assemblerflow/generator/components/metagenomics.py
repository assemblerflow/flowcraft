
try:
    from generator.process import Process
except ImportError:
    from assemblerflow.generator.process import Process

class bowtie_host(Process):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.params = {
            "refIndex": {
                "default": "'/index_hg19/hg19'",
                "description": "Specifies the reference indexes to be provided "
                               "to bowtie2."
            }
        }

        self.secondary_inputs = [
            {
                "params": "refIndex",
                "channel": "IN_index_files = Channel.value(params.refIndex)"
            }
        ]

        self.directives = {
            "bowtie_host": {
                "container": "cimendes/bowtie2_hg19",
                "version": "2.0"
            }
        }

        self.status_channels = [
            "bowtie_host"
        ]

