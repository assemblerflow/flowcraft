# Changelog

## Upcoming release (`dev` branch)

### New components

- `maxbin2`: An automatic tool for binning metagenomic sequences.
- `bowtie2`: Align short paired-end sequencing reads to long reference
sequences.

### New recipes

- `plasmids`: A recipe to perform mapping, mash screen on reads
and also mash dist for assembly based approaches (all to detect
plasmdis). This also includes annotation with abricate for the assembly.
- `plasmids_mapping`: A recipe to perfmorm mapping for plasmids.
- `plasmids_mash`: A recipe to perform mash screen for plasmids.
- `plasmids_assembly`: A recipe to perform mash dist for plasmid
assemblies.

### Minor/Other changes

- Added "smart" check when the user provides a typo in pipeline string
for a given process, outputting some "educated" guesses to the the
terminal.
- Added "-cr" option to show current recipe `pipeline_string`.
- Changed the way recipes were being parsed by `proc_collector` for the
usage of `-l` and `-L` options.
- Added check for non-ascii characteres in colored_print.

### Bug fixes

- Fixed pipeline names that contain new line characters.
- **Template: sistr.nf**: Fixed comparison that determined process status.

## 1.2.0

### New components

- `card_rgi`: Anti-microbial gene screening for assemblies
- `filter_poly`: Runs PrinSeq on paired-end FastQ files to remove low complexity sequences
- `kraken`: Taxonomical identification of FastQ files
- `megahit`: Metagenomic assembler for paired-end FastQ files
- `metaprob`: Preforms read binning on metagenomic FastQ files
- `metamlst`: Checks the Sequence Type of metagenomic FastQ reads using Multilocus Sequence Typing
- `metaspades`: Metagenomic assembler for paired-end FastQ files
- `midas_species`: Taxonomical identification of FastQ files on the species level
- `remove host`: Read mapping with Bowtie2 against the host genome (default hg19) and removes the mapping reads
- `sistr`: Salmonella *in silico* typing component for assemblies. 

### Features

- Added `inspect` run mode to flowcraft for displaying the progress overview
  during a nextflow run. This run mode has an `overview` and `broadcast` options
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
