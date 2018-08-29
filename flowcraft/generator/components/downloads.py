
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
            }
        }

        self.directives = {"reads_download": {
            "cpus": 1,
            "memory": "'1GB'",
            "container": "flowcraft/getseqena",
            "version": "0.4.0-1"
        }}


class FasterqDump(Process):
    """Process template for fasterq-dump

    This process is set with:

        - ``input_type``: accessions
        - ``output_type`` fastq

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "accessions"
        self.output_type = "fastq"

        self.params = {
            "option_file": {
                "default": "false",
                "description": "Read more options and parameters from the file."
                           "Use to provide parameters to fasterq-dump"
            },
            "compress_fastq": {
                "default": "true",
                "description": "This option allow the users to define if they"
                               "want to compress the downloaded fastq files, "
                               "saving disk space. Default behavior is set"
                               "to compress the fastq files. If the user wants"
                               "to change this, set the variable to 'no'"
            }
        }

        self.directives = {"fasterqDump": {
            "cpus": 1,
            "memory": "'1GB'",
            "container": "flowcraft/sra-tools",
            "version": "2.9.1-1"
        }}

        self.status_channels = [
            "fasterqDump"
        ]
