
try:
    from generator.process import Process
except ImportError:
    from flowcraft.generator.process import Process


class IntegrityCoverage(Process):
    """Process template interface for first integrity_coverage process

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains two **secondary channel link starts**:

        - ``SIDE_phred``: Phred score of the FastQ files
        - ``SIDE_max_len``: Maximum read length
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.params = {
            "genomeSize": {
                "default": 1,
                "description":
                    "Genome size estimate for the samples in Mb. It is used to "
                    "estimate the coverage and other assembly parameters and"
                    "checks"
            },
            "minCoverage": {
                "default": 0,
                "description":
                    "Minimum coverage for a sample to proceed. By default it's set"
                    "to 0 to allow any coverage"
            }
        }

        self.link_start.extend(["SIDE_phred", "SIDE_max_len"])


class CheckCoverage(Process):
    """Process template interface for additional integrity_coverage process

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains one **secondary channel link start**:

        - ``SIDE_max_len``: Maximum read length

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.params = {
            "genomeSize": {
                "default": 2.1,
                "description":
                    "Genome size estimate for the samples. It is used to "
                    "estimate the coverage and other assembly parameters and"
                    "checks"
            },
            "minCoverage": {
                "default": 15,
                "description":
                    "Minimum coverage for a sample to proceed. Can be set to"
                    "0 to allow any coverage"
            }
        }

        self.link_start.extend(["SIDE_max_len"])


class TrueCoverage(Process):
    """TrueCoverage process template interface
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.params = {
            "species": {
                "default": "null",
                "description":
                    "Species name. Must be the complete species name with"
                    "genus and species, e.g.: 'Yersinia enterocolitica'. "
            }
        }

        self.directives = {
            "true_coverage": {
                "cpus": 4,
                "memory": "'1GB'",
                "container": "flowcraft/true_coverage",
                "version": "3.2-1"
            }
        }


class FastQC(Process):
    """FastQC process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains two **status channels**:

        - ``STATUS_fastqc``: Status for the fastqc process
        - ``STATUS_report``: Status for the fastqc_report process

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.status_channels = ["STATUS_fastqc2", "STATUS_fastqc2_report"]
        """
        list: Setting status channels for FastQC execution and FastQC report
        """

        self.params = {
            "adapters": {
                "default": "'None'",
                "description":
                    "Path to adapters files, if any."
            }
        }

        self.directives = {"fastqc2": {
            "cpus": 2,
            "memory": "'4GB'",
            "container": "flowcraft/fastqc",
            "version": "0.11.7-1"
        }}


class Trimmomatic(Process):
    """Trimmomatic process template interface

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains one **secondary channel link end**:

        - ``SIDE_phred`` (alias: ``SIDE_phred``): Receives FastQ phred score
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_end.append({"link": "SIDE_phred", "alias": "SIDE_phred"})

        self.dependencies = ["integrity_coverage"]

        self.params = {
            "adapters": {
                "default": "'None'",
                "description":
                    "Path to adapters files, if any."
            },
            "trimSlidingWindow": {
                "default": "'5:20'",
                "description":
                    "Perform sliding window trimming, cutting once the "
                    "average quality within the window falls below a "
                    "threshold"
            },
            "trimLeading": {
                "default": "3",
                "description":
                    "Cut bases off the start of a read, if below a threshold "
                    "quality"
            },
            "trimTrailing": {
                "default": "3",
                "description":
                    "Cut bases of the end of a read, if below a "
                    "threshold quality"
            },
            "trimMinLength": {
                "default": "55",
                "description":
                    "Drop the read if it is below a specified length "
            }
        }

        self.directives = {"trimmomatic": {
            "cpus": 2,
            "memory": "{ 4.GB * task.attempt }",
            "container": "flowcraft/trimmomatic",
            "version": "0.36-1"
        }}


class FastqcTrimmomatic(Process):
    """Fastqc + Trimmomatic process template interface

    This process executes FastQC only to inform the trim range for trimmomatic,
    not for QC checks.

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    It contains one **secondary channel link end**:

        - ``SIDE_phred`` (alias: ``SIDE_phred``): Receives FastQ phred score

    It contains three **status channels**:

        - ``STATUS_fastqc``: Status for the fastqc process
        - ``STATUS_report``: Status for the fastqc_report process
        - ``STATUS_trim``: Status for the trimmomatic process
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_end.append({"link": "SIDE_phred", "alias": "SIDE_phred"})

        self.status_channels = ["STATUS_fastqc", "STATUS_fastqc_report",
                                "STATUS_trimmomatic"]

        self.dependencies = ["integrity_coverage"]

        self.params = {
            "adapters": {
                "default": "'None'",
                "description":
                    "Path to adapters files, if any."
            },
            "trimSlidingWindow": {
                "default": "'5:20'",
                "description":
                    "Perform sliding window trimming, cutting once the "
                    "average quality within the window falls below a "
                    "threshold."
            },
            "trimLeading": {
                "default": "3",
                "description":
                    "Cut bases off the start of a read, if below a threshold "
                    "quality."
            },
            "trimTrailing": {
                "default": "3",
                "description":
                    "Cut bases of the end of a read, if below a "
                    "threshold quality."
            },
            "trimMinLength": {
                "default": "55",
                "description":
                    "Drop the read if it is below a specified length."
            }
        }

        self.directives = {
            "fastqc": {
                "cpus": 2,
                "memory": "'4GB'",
                "container": "flowcraft/fastqc",
                "version": "0.11.7-1"
            },
            "trimmomatic": {
                "cpus": 2,
                "memory": "{ 4.GB * task.attempt }",
                "container": "flowcraft/trimmomatic",
                "version": "0.36-1"
            }
        }


class FilterPoly(Process):
    """PrinSeq process to filter non-informative sequences from reads

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq
        - ``ptype``: pre_assembly

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.params = {
            "adapter": {
                "default": "'A 50%; T 50%; N 50%'",
                "description":
                    "Pattern to filter the reads. Please separate parameter"
                    "values with a space and separate new parameter sets with semicolon (;)."
                    "Parameters are defined by two values: the pattern (any combination of the"
                    "letters ATCGN), and the number of repeats or percentage of occurence."
            }
        }

        self.directives = {"filter_poly": {
            "cpus": 1,
            "memory": "{ 4.GB * task.attempt }",
            "container": "flowcraft/prinseq",
            "version": "0.20.4-1"
        }}

        self.status_channels = [
            "filter_poly"
        ]


class DownsampleFastq(Process):
    """Downsamples FastQ file based on depth using seqtk

    This process is set with:

        - ``input_type``: fastq
        - ``output_type``: fastq

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.params = {
            "genomeSize": {
                "default": 1,
                "description":
                    "Genome size estimate for the samples in Mb. It is used to"
                    " estimate the coverage"
            },
            "depth": {
                "default": 100,
                "description":
                    "Maximum estimated depth coverage allowed. FastQ with "
                    "higher estimated depth will be subsampled to this value."
            }
        }

        self.directives = {"downsample_fastq": {
            "cpus": 1,
            "memory": "{ 4.GB * task.attempt }",
            "container": "flowcraft/seqtk",
            "version": "1.3.0-3"
        }}

        self.status_channels = [
            "downsample_fastq"
        ]