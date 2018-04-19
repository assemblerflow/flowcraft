try:
    from generator.process import Process
except ImportError:
    from assemblerflow.generator.process import Process


class PatlasMapping(Process):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.input_type = "fastq"
        self.output_type = "json"

        self.params = {
            "max_k": {
                "default": 10949,
                "description": "Sets the k parameter for bowtie2 allowing to "
                               "make multiple mappings of the same read "
                               "against several hits on the query sequence or "
                               "sequences."
            },
            "trim5": {
                "default": 0,
                "description": "Sets trim5 option for bowtie. This will become"
                               " legacy with QC integration, but it enables to"
                               " trim 5' end of reads to be mapped with "
                               "bowtie2."
            },
            "cov_cutoff": {
                "default": 0.6,
                "description": "This variable sets a cutoff for the percentage"
                               " of the query reference sequence that is "
                               "covered by reads (in absolute lenght)."
            },
            "plasmid_length_dict": {
                "default": "'jsons/*_length.json'",
                "description": "Stores a dictionary of lengths to be added to "
                               "jsonDumpingMapping process so that it can "
                               "easily get the size of each plasmid sequence "
                               "to be queried."
            },
            "refIndex": {
                "default": "'/ngstools/data/indexes/bowtie2idx/bowtie2.idx'",
                "description": "Specifies the reference indexes to be provided"
                               " to bowtie2."
            },
            "samtoolsIndex": {
                "default": "'/ngstools/data/indexes/fasta/samtools.fasta.fai'",
                "description": "Specifies the reference indexes to be provided"
                               " to samtools."
            },
            "lengthJson": {
                "default": "'/ngstools/data/reads_sample_result_length.json'",
                "description": "A dictionary of all the lengths of reference "
                               "sequences."
            }
        }

        self.secondary_inputs = [
            {
                "params": "refIndex",
                "channel": "IN_index_files = Channel.value(params.refIndex)"
            },
            {
                "params": "samtoolsIndex",
                "channel": "IN_samtools_indexes = Channel"
                           ".value(params.samtoolsIndex)"
            },
            {
                "params": "lengthJson",
                "channel": "IN_length_json = Channel.value(params.lengthJson)"
            }
        ]

        self.directives = {
            "mappingBowtie": {
                "container": "tiagofilipe12/patlasflow_mapping",
                "version": "1.1.2"
            },
            "samtoolsView": {
                "container": "tiagofilipe12/patlasflow_mapping",
                "version": "1.1.2"
            },
            "jsonDumpingMapping": {
                "container": "tiagofilipe12/patlasflow_mapping",
                "version": "1.1.2"
            }
        }

        self.status_channels = [
            "mappingBowtie",
            "samtoolsView",
            "jsonDumpingMapping"
        ]

        self.compiler["patlas_consensus"] = ["mappingOutputChannel"]
