# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

Snakemake workflow for RNA-Seq alignment using STAR, with duplicate marking (Picard/umi_tools), coverage track generation (deepTools), and MultiQC reporting. Integrates with BioRoot utilities (`bioroots_utilities.smk`) for sample management and organism configuration.

## Key Commands

```bash
# Run workflow with N cores
snakemake --cores <N>

# Run workflow with conda environments
snakemake --cores <N> --use-conda

# Run specific target
snakemake --cores <N> qc_reports/all_samples/alignment_RNA_multiqc/multiqc.html

# Dry run to check workflow
snakemake --cores <N> --dryrun
```

## Architecture

### Main Files
- **[Snakefile](e:\github\BioIT-CEITEC\alignment_rna\Snakefile)** - Entry point; loads config, imports BioRoot modules, defines `rule all` targets
- **[workflow.config.json](e:\github\BioIT-CEITEC\alignment_rna\workflow.config.json)** - Workflow metadata, GUI parameters, and user-configurable options
- **[rules/alignment_RNA.smk](e:\github\BioIT-CEITEC\alignment_rna\rules\alignment_RNA.smk)** - Core rules: `alignment_RNA`, `mark_duplicates`, `get_cov_tracks`, `alignment_RNA_multiqc`

### Wrapper Scripts (in `wrappers/`)
Each rule has a corresponding Python wrapper script that executes the actual tools:
- **alignment_RNA/script.py** - STAR alignment, BAM sorting, transcriptome BAM, bedGraph/bigWig generation, chimeric read handling
- **mark_duplicates/script.py** - Picard MarkDuplicates or umi_tools dedup depending on UMI settings
- **get_cov_tracks/script.py** - bamCoverage from deepTools for bedGraph/bigWig
- **alignment_RNA_multiqc/script.py** - MultiQC report aggregation

### Conda Environments
Each wrapper has its own environment defined in `env.yaml` files under their respective directories.

## Configuration Flow

1. **BioRoot sample loading**: `BR.load_sample()` reads sample metadata
2. **Organism config**: `BR.load_organism()` loads organism-specific paths (FASTA, GTF, STAR index, etc.)
3. **Read pair handling**: `BR.set_read_pair_qc_tags()`, `BR.set_read_pair_tags()`, `BR.set_paired_tags()` determine SE vs PE handling
4. **User parameters**: From [workflow.config.json](file://e:\github\BioIT-CEITEC\alignment_rna\workflow.config.json) - strandness, mismatch rates, UMI settings, duplicate removal, etc.

## Key Parameters (from workflow.config.json)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `RNAseq_type` | Library prep type (classic_fwd/rev, quant_fwd/rev, sense, corall) | classic_rev |
| `strandness` | unstr/fwd/rev | rev |
| `mark_duplicates` | Enable duplicate marking | true |
| `umi_usage` | UMI handling mode | mark_duplicates |
| `umi_sep` | UMI separator character | _ |
| `remove_duplicates` | Actually remove vs just mark | false |
| `map_score` / `map_perc` | Mapping quality thresholds | 0.66 |
| `perc_mismatch` | Max mismatch ratio | 0.1 |
| `max_intron` | Max intron size | 1000000 |
| `get_cov_tracks` | Generate coverage tracks | false |
| `keep_not_markDups_bam` | Keep pre-dedup BAM | false |

## Output Structure

```
mapped/
  {sample}.bam                    # Final aligned BAM
  {sample}.bam.bai                # BAM index
  {sample}/                       # STAR intermediate files
  transcriptome/{sample}.transcriptome.bam
  {sample}.bigWig                 # Coverage (if enabled)
  {sample}.bedGraph               # Coverage (if enabled)
qc_reports/
  {sample}/MarkDuplicates/{sample}.markDups_metrics.txt
  all_samples/alignment_RNA_multiqc/multiqc.html
logs/
  {sample}/alignment.log
  {sample}/mark_duplicates.log
  {sample}/get_cov_tracks.log
```

## Important Implementation Details

- STAR runs in **two-pass mode** (`--twopassMode Basic`) for improved junction detection
- Chimeric reads are converted to BAM separately from the main alignment
- Chromosome names are converted to UCSC format for bigWig visualization using [convert_chromosome_names.R](file://e:\github\BioIT-CEITEC\alignment_rna\wrappers\alignment_RNA\convert_chromosome_names.R)
- When UMI deduplication is enabled, `umi_tools dedup` with directional method replaces Picard
- Coverage tracks use bin size 5bp (`-bs 5`) via bamCoverage
- Memory resources are specified per-rule (alignment: 34GB, mark_duplicates: 15GB)
- Thread allocation: alignment (40), mark_duplicates (8), get_cov_tracks (4)
