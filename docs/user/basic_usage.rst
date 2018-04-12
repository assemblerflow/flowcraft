Basic Usage
===========

Assemblerflow has currently one execution mode, ``build``, that is used to
build the nextflow pipeline. However, more execution modes are slated for
release.

Assembling a pipeline
---------------------

Pipelines can be generated using the ``build`` execution mode of assemblerflow
and the ``-t`` parameter to specify the components inside quotes::

    assemblerflow build -t "trimmomatic fastqc spades" -o my_pipe.nf

This command will generate a linear pipeline with three components on the
current working directory (for more features and tips on how pipelines can be
built, see the :doc:`pipeline_building` section). A linear pipeline means that
there are no bifurcations between components, and the input data will flow
linearly across components. In this particular case, the input data of the
pipeline will be paired-end FastQ data, since that is the input data type
of the first component, ``trimmomatic``.

The rational of how the data flows across the pipeline is simple and intuitive.
Data enters a component and is processed in some way, which may result on the
creation of results (stored in the ``results`` directory) and reports (stored
in the ``reports`` directory). If that component has an ``output_type``, it
will feed the processed data into the next component (or components). This
will

Pipeline directory
------------------

In addition to the main nextflow pipeline file (``my_pipe.nf``),
assemblerflow will write several auxiliary files that are necessary for
the pipeline to run. The contents of the directory should look something like
this::

    $ ls
    bin                lib           my_pipe.nf       params.config     templates
    containers.config  my_pipe.html  nextflow.config  resources.config  user.config



Parameters
----------

The parameters of the pipeline can be viewed by running the pipeline file
with ``nextflow`` and using the ``--help`` option::

    $ nextflow my_pipe.nf --help

    N E X T F L O W  ~  version 0.28.0
    Launching `my_pipe.nf` [stupefied_booth] - revision: 504208431f

    ============================================================
                    A S S E M B L E R F L O W
    ============================================================
    Built using assemblerflow v1.0.2


    Usage:
        nextflow run my_pipe.nf

           --fastq                     Path expression to paired-end fastq files. (default: fastq/*_{1,2}.*) (integrity_coverage)
           --genomeSize                Genome size estimate for the samples. It is used to estimate the coverage and other assembly parameters andchecks (default: 2.1) (integrity_coverage)
           --minCoverage               Minimum coverage for a sample to proceed. Can be set to0 to allow any coverage (default: 15) (integrity_coverage)
           --adapters                  Path to adapters files, if any (default: None) (trimmomatic;fastqc)
           --trimSlidingWindow         Perform sliding window trimming, cutting once the average quality within the window falls below a threshold (default: 5:20) (trimmomatic)
           --trimLeading               Cut bases off the start of a read, if below a threshold quality (default: 3 (trimmomatic)
           --trimTrailing              Cut bases of the end of a read, if below a threshold quality (default: 3) (trimmomatic)
           --trimMinLength             Drop the read if it is below a specified length (default: 55) (trimmomatic)
           --spadesMinCoverage         The minimum number of reads to consider an edge in the de Bruijn graph during the assembly (default: 2) (spades)
           --spadesMinKmerCoverage     Minimum contigs K-mer coverage. After assembly only keep contigs with reported k-mer coverage equal or above this value (default: 2) (spades)
           --spadesKmers               If 'auto' the SPAdes k-mer lengths will be determined from the maximum read length of each assembly. If 'default', SPAdes will use the default k-mer lengths. (default: auto) (spades)

All these parameters are related to the components of the pipeline. However,
the main input parameter (or parameters) of the pipeline are always available.
Since this pipeline started with FastQ paired-end files as the main input,
the ``--fastq`` parameter is available. If the pipeline started with any other
input type or with more than one input type, the appropriate parameters would
become available.

Executing the pipeline
----------------------

Most parameters in assemblerflow's components already come with sensible
defaults, which means that usually you'll only need to provide a small number
of arguments. In the example above, the ``--fastq`` is the only parameter
required. I have placed fastq files on the ``data`` directory, so::

    nextflow run my_pipe.nf --fastq "data/*_{1,2}.*"

The pattern for the fastq files is perhaps a bit confusing at first, but it's
necessary for the correct inference of the pairs. My fastq files are::

    $ ls data
    sample_1.fastq.gz  sample_2.fastq.gz

In this case, the pattern is given by the "_1." or "_2." substring, which leads
to the pattern ``*_{1,2}.*``. Another common nomenclature for paired fastq
files is something like ``sample_R1_L001.fastq.gz``. In this case, an
acceptable pattern would be ``*_R{1,2}_*``.

.. important::

    Note the quotes around the fastq path pattern. These quotes are necessary
    to allow nextflow to resolve the pattern, otherwise your shell might try
    to resolve the pattern and provide the wrong input to nextflow.

Results
-------


