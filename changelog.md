# Changelog

## Upcoming release (`dev` branch)

### New components

- `maxbin2`: An automatic tool for binning metagenomic sequences
- `bowtie2`: Align short paired-end sequencing reads to long reference sequences

### Minor/Other changes

- Added "smart" check when the user provides a typo in pipeline string
for a given process, outputting some "educated" guesses to the the
terminal.

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
