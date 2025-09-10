# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a single-cell RNA sequencing (scRNA-seq) analysis pipeline built on Seurat v5 for R>=4.2. The pipeline processes 10X Genomics Cell Ranger output and supports both regular samples and xenograft samples (human samples with mouse contamination). The workflow includes quality control, filtering, integration/merging, clustering, and downstream analysis.

## Key Commands

### Main Pipeline Execution
- `./CMDS.p12 [options] <parameters>` - Main pipeline script
  - `-i|--integrate` - Use integration instead of simple merge for combining samples
  - `-c|--config config.yaml` - Use explicit config file for multiplex runs
  - `-f|--filter` - Filter only (stop after stage 2a)
- `./getCRConfigTmpl.sh <CELL_RANGER_DIR>` - Generate config template for multiplex runs

### Pipeline Stages
The pipeline runs in sequential stages:
1. **Stage 1**: `doSeuratV5_01.R` - Initial QC and filtering parameter determination
2. **Stage 2a**: `doSeuratV5_02a_MergeOnly.R` or `doSeuratV5_02a_IntegrateData.R` - Sample combination
3. **Stage 2b**: `doSeuratV5_02b.R` - PCA and clustering
4. **Stage 3**: Various downstream analysis scripts (`doSeuratV5_03_*.R`)

## Architecture and Workflow

### Input Data Handling
- **Standard runs**: Auto-discovery from `../cellranger/s_*` directories
- **Multiplex runs**: Requires explicit `config.yaml` with sample paths and IDs
- **Xenograft support**: Uses `refdata-gex-GRCh38-and-mm10-2020-A` genome for human/mouse mixed samples

### Configuration System
- `pass_00_PARAMS.yaml` - Parameter overrides (editable)
- `data/genomes.yaml` - Genome mapping definitions
- Generated parameter files: `pass_01_PARAMS.yaml`, `pass_02_PARAMS.yaml`

### Core R Modules
- `seuratTools.R` - Main Seurat utilities and helper functions
- `qcAndFilter.R` - Quality control and filtering functions
- `plotTools.R` - Plotting utilities
- `rsrc/` - Additional plotting and QC resources

### Analysis Features
- Quality control with MAD-based filtering (`MED +/- 3*MAD` criteria)
- Cell cycle regression analysis and plots
- Gene filtering in Stage 2a
- Support for both merge and integration workflows
- Downstream analysis: cell type annotation, feature plots, heatmaps, module scoring

### Project Structure
- Main analysis scripts: `doSeuratV5_*.R`
- Utilities: `*Tools.R`, `filterGenes.R`
- Resources: `rsrc/`, `opt/`, `data/`
- Results tracking: `RESULTS.md`, `RESULTS.pdf`
- Configuration: `config_example.yaml`, parameter YAML files

## Development Notes

- Project naming: Auto-extracted from path using `extractProjectIDFromPath.py`
- Debugging: Use `DEBUG=TRUE` for 10% downsampling during development
- The pipeline generates parameter files between stages that carry forward settings
- Cell cycle genes are defined in `data/cc.genes.mouse.v2.yaml`
- Supports both mouse (mm10) and human (hg38) genomes plus xenograft analysis