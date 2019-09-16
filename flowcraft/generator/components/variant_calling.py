try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class Haplotypecaller(Process):
    """Call germline SNPs and indels via local re-assembly of haplotypes

        This process is set with:

            - ``input_type``: bam
            - ``output_type``: vcf
            - ``ptype``: varaint calling

        """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "bam"

        self.params = {
            "reference": {
                "default": "null",
                "description": "Specifies the reference genome to be provided "
                               "to GATK HaplotypeCaller."
            },
            "intervals": {
                "default": "null",
                "description": "Interval list file to specify the regions to call variants."
            }
        }

        self.directives = {
            "haplotypecaller": {
                "container": "broadinstitute/gatk",
                "memory": "{2.Gb*task.attempt}",
                "cpus": 4,
            },
            "merge_vcfs": {
                "container": "broadinstitute/gatk",
                "memory": "{5.Gb*task.attempt}",
                "cpus": 4,
            }
        }

        self.status_channels = [
            "haplotypecaller",
            "merge_vcfs"
        ]
