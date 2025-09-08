# RNA-Seq Alignment Workflow

This repository contains a Snakemake workflow for processing and aligning RNA-Seq data. The workflow is modular, configurable, and designed for reproducible analysis in bioinformatics projects.

## Features
- **Modular Design:** Each step is implemented as a separate rule and wrapper, allowing easy customization.
- **Conda Environments:** Each rule uses its own environment for reproducibility.
- **BioRoot Utilities:** Integrates with external modules for sample and reference management.
- **Configurable:** All parameters (paths, organism, alignment options, etc.) are set via `workflow.config.json`.

## Workflow Overview
1. **Preprocessing** (optional):
   - Trimming and cleaning of raw FASTQ files.
   - Adapter removal and quality filtering.
2. **Quality Control:**
   - FastQC reports for cleaned FASTQ files.
   - MultiQC summary report for all samples.
3. **Alignment:**
   - STAR alignment of reads to reference genome.
   - Generation of BAM files and transcriptome BAM files.
4. **Mark Duplicates:**
   - Identification and marking/removal of duplicate reads.
   - Metrics reporting.
5. **Coverage Tracks:**
   - Generation of BigWig and BedGraph coverage files from BAMs.

## Directory Structure
```
Snakefile
workflow.config.json
rules/
    alignment_RNA.smk
wrappers/
    alignment_RNA/
        convert_chromosome_names.R
        env.yaml
        script.py
    alignment_RNA_multiqc/
        env.yaml
        multiqc_config.txt
        script.py
    cleaned_fastq_qc/
        env.yaml
        script.py
    get_cov_tracks/
        env.yaml
        script.py
    mark_duplicates/
        env.yaml
        script.py
    preprocess/
        env.yaml
        script.py
```

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
