# Changelog

## Upcoming release (`dev` branch)

## 1.2.1

### Improvements

- The parameter system has been revamped, and parameters are now component-specific
and independent by default. This allows a better fine-tuning of the parameters
and also the execution of the same component multiple times (for instance in a fork)
with different parameters. The old parameter system that merged identical parameters
is still available by using the `--merge-params` flag when building the pipeline.
- Added a global `--clearAtCheckpoint` parameter that, when set to true, will remove
temporary files that are no longer necessary for downstream steps of the pipeline
from the work directory. This option is currently supported for the `trimmomatic`,
`fastqc_trimmomatic`, `skesa` and `spades` components. 

### New components

- `maxbin2`: An automatic tool for binning metagenomic sequences.
- `bowtie2`: Align short paired-end sequencing reads to long reference
sequences.
- `retrieve_mapped`: Retrieves the mapped reads of a previous bowtie2 mapping process.

### New recipes

- `plasmids`: A recipe to perform mapping, mash screen on reads
and also mash dist for assembly based approaches (all to detect
plasmids). This also includes annotation with abricate for the assembly.
- `plasmids_mapping`: A recipe to perform mapping for plasmids.
- `plasmids_mash`: A recipe to perform mash screen for plasmids.
- `plasmids_assembly`: A recipe to perform mash dist for plasmid
assemblies.

### Minor/Other changes

- Added "smart" check when the user provides a typo in pipeline string
for a given process, outputting some "educated" guesses to the
terminal.
- Added "-cr" option to show current recipe `pipeline_string`.
- Changed the way recipes were being parsed by `proc_collector` for the
usage of `-l` and `-L` options.
- Added check for non-ascii characters in colored_print.
- Fixed log when a file with the pipeline is provided to -t option
instead of a string.

### Bug fixes

- Fixed pipeline names that contain new line characters.
- Fixed pipeline generation when automatic dependencies were added right after a fork
- **Template: sistr.nf**: Fixed comparison that determined process status.
- Fixed issue with `--version` option.

## 1.2.0

### New components

- `card_rgi`: Anti-microbial resistance gene screening for assemblies
- `filter_poly`: Runs PrinSeq on paired-end FastQ files to remove low complexity sequences
- `kraken`: Taxonomic identification on FastQ files
- `megahit`: Metagenomic assembler for paired-end FastQ files
- `metaprob`: Performs read binning on metagenomic FastQ files
- `metamlst`: Checks the Sequence Type of metagenomic FastQ reads using Multilocus Sequence Typing
- `metaspades`: Metagenomic assembler for paired-end FastQ files
- `midas_species`: Taxonomic identification on FastQ files at the species level
- `remove host`: Read mapping with Bowtie2 against the target host genome (default hg19) and removes the mapping reads
- `sistr`: Salmonella *in silico* typing component for assemblies. 

### Features

- Added `inspect` run mode to flowcraft for displaying the progress overview
  during a nextflow run. This run mode has `overview` and `broadcast` options
  for viewing the progress of a pipeline.

### Minor/Other changes

- Changed `mapping_patlas` docker container tag and variable
(PR [#76](https://github.com/assemblerflow/assemblerflow/pull/76)).
- The `env` scope of nextflow.config now extends the `PYTHONPATH`
environmental variable.
- Updated indexes for both `mapping_patlas` and `mash` based processes.
- New logo!

### Bug Fixes

- **Template: fastqc_report.py**: Added fix to trim range evaluation.
- **Script: merge_json.py**: Fixed chewbbaca JSON merge function.
