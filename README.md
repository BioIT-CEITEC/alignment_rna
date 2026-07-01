# RNA-Seq Alignment Workflow

This repository contains a Snakemake workflow for processing and aligning RNA-Seq data. The workflow is modular, configurable, and designed for reproducible analysis in bioinformatics projects.

## Features
- **Modular Design:** Each step is implemented as a separate rule and wrapper, allowing easy customization.
- **Conda Environments:** Each rule uses its own environment for reproducibility.
- **BioRoot Utilities:** Integrates with external modules for sample and reference management.
- **Configurable:** All parameters (paths, organism, alignment options, etc.) are set via `workflow.config.json`.

## Workflow Overview
1. **Alignment:**
   - STAR alignment of reads to reference genome.
   - Generation of BAM files and transcriptome BAM files.
2. **Mark Duplicates:**
   - Identification and marking/removal of duplicate reads.
   - Metrics reporting.
3. **Coverage Tracks:**
   - Generation of BigWig and BedGraph coverage files from BAMs.

## Directory Structure
- `Snakefile`: Main workflow file.
- `workflow.config.json`: Configuration file.
- `rules/`: Snakemake rule files for each workflow step.
- `wrappers/`: Scripts and conda environments for each step.
- `qc_reports/`: Output directory for QC results and reports.
- `processed_fastq/`: Input directory for processed FastQ files.
- `mapped/`: Output directory for alignment BAM files.
- `logs/`: Log files for each step.

## Usage
1. **Configure the workflow:**
   - Edit `workflow.config.json` to specify sample information, reference genome, and parameters.
2. **Run the workflow:**
   ```bash
   snakemake --cores <N>
   ```
   Replace `<N>` with the number of CPU cores to use.
3. **Outputs:**
   - Aligned BAM files, coverage tracks, QC reports, and summary HTML files in the respective output directories.

## Requirements
- [Snakemake >=5.18.0](https://snakemake.readthedocs.io/)
- [Conda](https://docs.conda.io/)
- Python 3

## Customization
- Modify rules or wrapper scripts to adapt to specific project needs.
- Add or remove steps by editing the `rules/` and `wrappers/` directories.

## Contact
For questions or contributions, please contact the BioIT-CEITEC team.
