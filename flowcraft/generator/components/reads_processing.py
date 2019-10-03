
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Fq2FaPaired(Process):
    """Fq2Fa_Paired process template interface
        This process uses bbtools reformat.sh to convert
        paired end fastq reads to paired end fasta reads
        (without quality information

        This process is set with:

            - ``input_type``: fastq
            - ``output_type``: fasta_se
            - ``ptype``: reformat

        """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta_pe"

        self.directives = {"fq2fa_paired": {
            "container": "pcerqueira/bbtools",
            "version": "38.44"
        }}


class Fq2FaSingle(Process):
    """Fq2Fa_Paired process template interface
        This process uses bbtools reformat.sh to convert
        paired end fastq reads to a single fasta reads file
        (without quality information

        This process is set with:

            - ``input_type``: fastq
            - ``output_type``: fasta_se
            - ``ptype``: reformat

        """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fasta_se"

        self.directives = {"fq2fa_single": {
            "container": "pcerqueira/bbtools",
            "version": "38.44"
        }}
