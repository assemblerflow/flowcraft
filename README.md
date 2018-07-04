# FlowCraft :whale2::package:

![Nextflow version](https://img.shields.io/badge/nextflow->0.27.0-brightgreen.svg)
![Python version](https://img.shields.io/badge/python-3.6-brightgreen.svg)
[![Build Status](https://travis-ci.org/assemblerflow/flowcraft.svg?branch=master)](https://travis-ci.org/assemblerflow/flowcraft)
[![codecov](https://codecov.io/gh/assemblerflow/flowcraft/branch/master/graph/badge.svg)](https://codecov.io/gh/assemblerflow/flowcraft)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/f518854f780b41a08ca2fb1c14e360f0)](https://www.codacy.com/app/o.diogosilva/assemblerflow?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ODiogoSilva/assemblerflow&amp;utm_campaign=Badge_Grade)
[![Documentation Status](https://readthedocs.org/projects/flowcraft/badge/?version=latest)](http://flowcraft.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/flowcraft.svg)](https://badge.fury.io/py/flowcraft)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/flowcraft/badges/version.svg)](https://anaconda.org/bioconda/flowcraft)

<p align="center">
  <img width="360" src="docs/resources/logo_large.png" alt="nextflow_logo"/>
</p>

A [Nextflow](https://www.nextflow.io/) pipeline assembler for genomics.
Pick your modules. Assemble them. Run the pipeline.

(Previously known as Assemblerflow)

## The premisse

#### Build a pipeline

What if building your own genomics pipeline would be as simple as:

```
flowcraft.py build -t "trimmomatic fastqc skesa pilon" -o my_pipeline.nf
```

Seems pretty simple right? What if we could run this pipeline with a single command on any linux machine or cluster by leveraging
the awesomeness of [nextflow](https://www.nextflow.io/) and [docker](https://www.docker.com/)/[singularity](http://singularity.lbl.gov/)
containers without having to install any of the pipeline dependencies?

#### Run the pipeline

```
nextflow run my_pipeline.nf --fastq path/to/fastq

N E X T F L O W  ~  version 0.30.1
Launching `my_pipeline.nf` [admiring_lamarck] - revision: 82cc9cd2ed

============================================================
                M Y   P I P E L I N E
============================================================
Built using flowcraft v1.2.1

 Input FastQ                 : 2
 Input samples               : 1
 Reports are found in        : ./reports
 Results are found in        : ./results
 Profile                     : standard

Starting pipeline at Tue Jun 12 19:38:26 WEST 2018

[warm up] executor > local
[7c/eb5f2f] Submitted process > integrity_coverage_1_1 (02AR0553)
(...)
[31/7d90a1] Submitted process > compile_pilon_report_1_6

Completed at: Tue Jun 12 19:58:32 WEST 2018
Duration    : 20m 6s
Success     : true
Exit status : 0
```

Congratulations! You just built and executed your own pipeline with
only two commands! :tada:

## Installation

FlowCraft is available as a bioconda package, which already brings
nextflow:

```
conda install flowcraft
```

#### Container engines

Pipelines built with FlowCraft require at least one container
engine to be installed, among `docker`, `singularity` or `shifter`.
If you already have any one of these installed, you're good to go.
If not, we recommend installing singularity, which also has a bioconda
package:

```
conda install singularity
```

## How to use it

The complete user guide of FlowCraft can be found on [readthedocs.org](http://flowcraft.readthedocs.io/en/latest/?badge=latest).
For a quick and dirty demonstration, see below.

### Quick guide

#### Building a pipeline

FlowCraft comes with a number of [ready-to-use components](http://flowcraft.readthedocs.io/en/latest/user/available_components.html) to build your
own pipeline. Following some basic rules, such as the output type of one process
must match the input type of the next process, assembling a pipeline is done
using the `build` mode and the `-t` option:

```
flowcraft build -t "trimmomatic spades abricate" -o my_pipeline.nf -n "assembly pipe"
```

This command will generate everything that is necessary to run the
pipeline automatically, but the main pipeline executable
file will be `my_pipeline.nf`. This file will contain a nextflow pipeline
for genome assembly starts with `trimmomatic` and finishes with anti-microbial
gene annotation using `abricate`.

#### Wait... what about the software parameters?

Each component in the pipeline has its own set of parameters that can be
modified before or when executing the pipeline. These parameters are
described in the documentation of each process and you can check the options
of your particular pipeline using the `help` option:

```
$ nextflow run my_pipeline.nf --help
N E X T F L O W  ~  version 0.30.1
Launching `my_pipeline.nf` [prickly_picasso] - revision: 2e1a226e6d

============================================================
                F L O W C R A F T
============================================================
Built using flowcraft v1.2.1


Usage: 
    nextflow run my_pipeline.nf

       --fastq                     Path expression to paired-end fastq files. (default: fastq/*_{1,2}.*) (default: 'fastq/*_{1,2}.*')
       
       Component 'INTEGRITY_COVERAGE_1_1'
       ----------------------------------
       --genomeSize_1_1            Genome size estimate for the samples in Mb. It is used to estimate the coverage and other assembly parameters andchecks (default: 1)
       --minCoverage_1_1           Minimum coverage for a sample to proceed. By default it's setto 0 to allow any coverage (default: 0)
       
       Component 'TRIMMOMATIC_1_2'
       ---------------------------
       --adapters_1_2              Path to adapters files, if any. (default: 'None')
       --trimSlidingWindow_1_2     Perform sliding window trimming, cutting once the average quality within the window falls below a threshold (default: '5:20')
       --trimLeading_1_2           Cut bases off the start of a read, if below a threshold quality (default: 3)
       --trimTrailing_1_2          Cut bases of the end of a read, if below a threshold quality (default: 3)
       --trimMinLength_1_2         Drop the read if it is below a specified length  (default: 55)
       
       Component 'FASTQC_1_3'
       ----------------------
       --adapters_1_3              Path to adapters files, if any. (default: 'None')
       
       Component 'ASSEMBLY_MAPPING_1_5'
       --------------------------------
       --minAssemblyCoverage_1_5   In auto, the default minimum coverage for each assembled contig is 1/3 of the assembly mean coverage or 10x, if the mean coverage is below 10x (default: 'auto')
       --AMaxContigs_1_5           A warning is issued if the number of contigs is overthis threshold. (default: 100)
       --genomeSize_1_5            Genome size estimate for the samples. It is used to check the ratio of contig number per genome MB (default: 2.1)
```

This help message is dynamically generated depending on the pipeline you build.
Since this pipeline starts with `trimmomatic`, which receives fastq files as input,
`--fastq` is the default parameter for providing paired-end fastq files.

#### Running a pipeline

Now that we have our nextflow pipeline built, we are ready to executed it by
providing input data. By default, FlowCraft pipelines will run locally and use
`singularity` to run the containers of each component. This can be
changed in multiple ways, but for convenience FlowCraft has already defined
profiles for most configurations of `executors` and `container` engines.

Running a pipeline locally with `singularity` can be done with:

```
# Pattern for paired-end fastq is '<sample>_1.fastq.gz <sample>_2.fastq.gz'
nextflow run my_pipeline --fastq "path/to/fastq/*_{1,2}.*"
```

If you want to run a pipeline in a cluster with SLURM and singularity, just use
the appropriate profile:

```
nextflow run my_pipeline --fastq "path/to/fastq/*_{1,2}.*" -profile slurm_sing
```

During the execution of the pipeline, the results and reports for each component
are continuously saved to the `results` and `reports` directory, respectively.

#### Inspecting a pipeline progress

Since version 1.2.0, it is possible to inspect the progress of a nextflow pipeline
using the `flowcraft inspect` mode. To check the progress in a terminal, simply
type:

```
flowcraft inspect
```

On the directory where the pipeline is running. Alternatively, you can view the progress
in FlowCraft's web service by using the ``broadcast`` option:

```
flowcraft inspect -m broadcast
```

<img src="https://github.com/assemblerflow/flowcraft-webapp/raw/master/flowcraft-webapp/frontend/resources/fc_short_demo.gif"/>

## Why not just write a Nextflow pipeline?

In many cases, building a static nextflow pipeline is sufficient for our goals.
However, when building our own pipelines, we often felt the need to add dynamism
to this process, particularly if we take into account how fast new tools arise
and existing ones change. Our biological goals also change over time and we
might need different pipelines to answer different questions. FlowCraft makes
this very easy, by having a set of pre-made and ready-to-use components that can
be freely assembled.

For instance, changing the assembly software in a genome assembly pipeline becomes
as easy as:

```
# Use spades
trimmomatic spades pilon
# Use skesa
trimmomatic skesa pilon
```

![example1](https://github.com/ODiogoSilva/flowcraft/raw/master/docs/resources/example_3.png)

If you are interested in having some sort of genome annotation, simply add those
components at the end, using a fork syntax:

```
# Run prokka and abricate at the end of the assembly
trimmomatic spades pilon (prokka | abricate)
```

![example2](https://github.com/ODiogoSilva/flowcraft/raw/master/docs/resources/example_1.png)

On the other hand, if you are interest in just perform allele calling for wgMLST,
simply add `chewbbaca`:

```
trimmomatic spades pilon chewbbaca
```

![example3](https://github.com/ODiogoSilva/flowcraft/raw/master/docs/resources/example_2.png)

Since nextflow handles parallelism of large sets of data so well, simple pipelines
of two components are also useful to build:

```
trimmomatic fastqc
```

As the number of existing components grow, so does your freedom to build pipelines.

## Roadmap

You can see what we're planning next on our [roadmap guide](https://github.com/assemblerflow/flowcraft/wiki/Roadmap).

## Developer guide

### Adding new components

Is there a missing component that you would like to see included? We would love
to expand! You could make a component request in our
[issue tracker](https://github.com/assemblerflow/flowcraft/issues).

If you want to be part of the team, you can contribute with the code as well. Each component
in FlowCraft can be independently added without having to worry about
the rest of the code base. You'll just need to have some knowledge of python
and nextflow. [Check the developer documentation for how-to guides](http://assemblerflow.readthedocs.io/en/latest/)
