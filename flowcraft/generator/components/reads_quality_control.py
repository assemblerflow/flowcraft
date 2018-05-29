
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
                    "checks (default: $params.genomeSize)"
            },
            "minCoverage": {
                "default": 0,
                "description":
                    "Minimum coverage for a sample to proceed. By default it's set"
                    "to 0 to allow any coverage (default: $params.minCoverage)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "genomeSize",
                "channel":
                    "IN_genome_size = Channel"
                    ".value(params.genomeSize)"
                    "map{it -> it.toString().isNumber() ?"
                    " it : exit(1, \"The genomeSize parameter must be a number"
                    "or a float. Provided value: '${params.genomeSize}'\")}"
            },
            {
                "params": "minCoverage",
                "channel":
                    "IN_min_coverage = Channel"
                    ".value(params.minCoverage)"
                    "map{it -> it.toString().isNumber() ?"
                    " it : exit(1, \"The minCoverage parameter must be a "
                    "number or a float. Provided value: "
                    "'${params.minCoverage}'\")}"
            }
        ]

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
                    "checks (default: $params.genomeSize)"
            },
            "minCoverage": {
                "default": 15,
                "description":
                    "Minimum coverage for a sample to proceed. Can be set to"
                    "0 to allow any coverage (default: $params.minCoverage)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "genomeSize",
                "channel":
                    "IN_genome_size = Channel"
                    ".value(params.genomeSize)"
                    "map{it -> it.toString().isNumber() ?"
                    " it : exit(1, \"The genomeSize parameter must be a number"
                    "or a float. Provided value: '${params.genomeSize}'\")}"
            },
            {
                "params": "minCoverage",
                "channel":
                    "IN_min_coverage = Channel"
                    ".value(params.minCoverage)"
                    "map{it -> it.toString().isNumber() ?"
                    " it : exit(1, \"The minCoverage parameter must be a "
                    "number or a float. Provided value: "
                    "'${params.minCoverage}'\")}"
            }
        ]

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
                    "(default: $params.species)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "species",
                "channel":
                    "if ( !params.species){ exit 1, \"'species' parameter "
                    "missing\" }\n"
                    "if ( params.species.toString().split(\" \").size() != 2 )"
                    "{ exit 1, \"'species' parameter must contain two "
                    "values (e.g.: 'escherichia coli').Provided value: "
                    "'${params.species}'\"}\n"
                    "IN_pathoSpecies = Channel.value(params.species)"
            }
        ]

        self.directives = {
            "true_coverage": {
                "cpus": 4,
                "memory": "'1GB'",
                "container": "odiogosilva/true_coverage",
                "version": "3.2"
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
                    "Path to adapters files, if any "
                    "(default: $params.adapters)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "adapters",
                "channel": "IN_adapters = Channel.value(params.adapters)"
            }
        ]

        self.directives = {"fastqc2": {
            "cpus": 2,
            "memory": "'4GB'",
            "container": "ummidock/fastqc",
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
                    "Path to adapters files, if any "
                    "(default: $params.adapters)"
            },
            "trimSlidingWindow": {
                "default": "'5:20'",
                "description":
                    "Perform sliding window trimming, cutting once the "
                    "average quality within the window falls below a "
                    "threshold (default: $params.trimSlidingWindow)"
            },
            "trimLeading": {
                "default": "3",
                "description":
                    "Cut bases off the start of a read, if below a threshold "
                    "quality (default: $params.trimLeading"
            },
            "trimTrailing": {
                "default": "3",
                "description":
                    "Cut bases of the end of a read, if below a "
                    "threshold quality (default: $params.trimTrailing)"
            },
            "trimMinLength": {
                "default": "55",
                "description":
                    "Drop the read if it is below a specified length "
                    "(default: $params.trimMinLength)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "trimOpts",
                "channel":
                    "// Check sliding window parameter\n"
                    "if ( params.trimSlidingWindow.toString().split(\":\")"
                    ".size() != 2 )"
                    "{ exit 1, \"'trimSlidingWindow' parameter must contain"
                    "two values separated by a ':'. Provided value: "
                    "'${params.trimSlidingWindow}' \"}\n"
                    "if ( !params.trimLeading.toString().isNumber() )"
                    "{ exit 1, \"'trimLeading' parameter must be a number."
                    "Provide value: '${params.trimLeading}'\"}\n"
                    "if ( !params.trimTrailing.toString().isNumber() )"
                    "{ exit 1, \"'trimTrailing' parameter must be a number."
                    "Provide value: '${params.trimTrailing}'\"}\n"
                    "if ( !params.trimMinLength.toString().isNumber() )"
                    "{ exit 1, \"'trimMinLength' parameter must be a number."
                    "Provide value: '${params.trimMinLength}'\"}\n"
                    "IN_trimmomatic_opts = Channel."
                    "value([params.trimSlidingWindow,"
                    "params.trimLeading,params.trimTrailing,"
                    "params.trimMinLength])"
            },
            {
                "params": "adapters",
                "channel": "IN_adapters = Channel.value(params.adapters)"
            }
        ]

        self.directives = {"trimmomatic": {
            "cpus": 2,
            "memory": "{ 4.GB * task.attempt }",
            "container": "ummidock/trimmomatic",
            "version": "0.36-2"
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
                    "Path to adapters files, if any "
                    "(default: $params.adapters)"
            },
            "trimSlidingWindow": {
                "default": "'5:20'",
                "description":
                    "Perform sliding window trimming, cutting once the "
                    "average quality within the window falls below a "
                    "threshold (default: $params.trimSlidingWindow)"
            },
            "trimLeading": {
                "default": "3",
                "description":
                    "Cut bases off the start of a read, if below a threshold "
                    "quality (default: $params.trimLeading"
            },
            "trimTrailing": {
                "default": "3",
                "description":
                    "Cut bases of the end of a read, if below a "
                    "threshold quality (default: $params.trimTrailing)"
            },
            "trimMinLength": {
                "default": "55",
                "description":
                    "Drop the read if it is below a specified length "
                    "(default: $params.trimMinLength)"
            }
        }

        self.secondary_inputs = [
            {
                "params": "adapters",
                "channel": "IN_adapters = Channel.value(params.adapters)"
            },
            {
                "params": "trimOpts",
                "channel":
                    "// Check sliding window parameter\n"
                    "if ( params.trimSlidingWindow.toString().split(\":\")"
                    ".size() != 2 )"
                    "{ exit 1, \"'trimSlidingWindow' parameter must contain"
                    "two values separated by a ':'. Provided value: "
                    "'${params.trimSlidingWindow}' \"}\n"
                    "if ( !params.trimLeading.toString().isNumber() )"
                    "{ exit 1, \"'trimLeading' parameter must be a number."
                    "Provide value: '${params.trimLeading}'\"}\n"
                    "if ( !params.trimTrailing.toString().isNumber() )"
                    "{ exit 1, \"'trimTrailing' parameter must be a number."
                    "Provide value: '${params.trimTrailing}'\"}\n"
                    "if ( !params.trimMinLength.toString().isNumber() )"
                    "{ exit 1, \"'trimMinLength' parameter must be a number."
                    "Provide value: '${params.trimMinLength}'\"}\n"
                    "IN_trimmomatic_opts = Channel."
                    "value([params.trimSlidingWindow,"
                    "params.trimLeading,params.trimTrailing,"
                    "params.trimMinLength])"
            }
        ]

        self.directives = {
            "fastqc": {
                "cpus": 2,
                "memory": "'4GB'",
                "container": "ummidock/fastqc",
                "version": "0.11.7-1"
            },
            "trimmomatic": {
                "cpus": 2,
                "memory": "{ 4.GB * task.attempt }",
                "container": "ummidock/trimmomatic",
                "version": "0.36-2"
            }
        }
