
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class DownloadReads(Process):
    """Process template interface for reads downloading from SRA and NCBI

    This process is set with:

        - ``input_type``: accessions
        - ``output_type`` fastq

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "accessions"
        self.output_type = "fastq"

        self.params = {
            "asperaKey": {
                "default": "null",
                "description":
                    "Downloads fastq accessions from ENA using Aspera Connect "
                    "by providing the private-key file "
                    "'asperaweb_id_dsa.openssh' normally found in "
                    "~/.aspera/connect/etc/asperaweb_id_dsa.openssh "
                    "(Default: null)"
            }
        }

        self.directives = {"reads_download": {
            "cpus": 1,
            "memory": "'1GB'",
            "container": "ummidock/getseqena",
            "version": "0.4.0-2"
        }}
