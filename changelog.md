# Changelog

## New in upcoming 1.2.0 release (`dev` branch)

### New components

- `sistr`: Salmonella *in silico* typing component for assemblies. 

### Features

- Added `inspect` run mode to assemblerflow for displaying the progress overview
  during a nextflow run.

### Minor/Other changes

- Changed `mapping_patlas` docker container tag and variable
(PR [#76](https://github.com/assemblerflow/assemblerflow/pull/76)).
- The `env` scope of nextflow.config now extends the `PYTHONPATH` environmental variable

# Bug Fixes

- **Template: fastqc_report.py**: Added fix to trim range evaluation.
- **Script: merge_json.py**: Fixed chewbbaca JSON merge function.