# Codebase Overview

ampliconseq is a Nextflow pipeline for calling SNVs and indels in targeted amplicon sequencing data, developed at the Cancer Research UK Cambridge Institute.

## Project Structure

```
ampliconseq.nf          # Main Nextflow pipeline (DSL 2)
nextflow.config          # Pipeline configuration and execution profiles
pom.xml                  # Maven build for Java tools (Java 17, htsjdk)
download_vep_cache.nf    # Supplementary workflow for VEP cache download

bin/                     # R scripts (~20) for data processing and reporting
src/                     # Java source for custom bioinformatics tools
templates/               # Shell script templates for Nextflow processes
shiny/                   # R/Shiny web app for variant visualization
docker/                  # Dockerfile and conda environment spec
resources/               # Default reference data (blacklists, specific variants)
```

## Pipeline Architecture

The main workflow in `ampliconseq.nf` runs these stages:

1. **Input validation** - check sample sheet and amplicon coordinates
2. **Amplicon grouping** - create non-overlapping amplicon groups for variant calling
3. **Alignment metrics** - Picard CollectAlignmentSummaryMetrics and CollectTargetedPcrMetrics
4. **Read extraction** - extract reads per amplicon group (Java tool)
5. **Variant calling** - VarDict (default), GATK HaplotypeCaller, or GATK Mutect2
6. **Pileup counts** - allele counts at each position (Java tool)
7. **Collation** - merge VCFs, convert to tabular format
8. **Background noise filtering** - fit Beta distributions to model substitution noise
9. **VEP annotation** (optional) - Ensembl Variant Effect Predictor
10. **Summary** - final variant table with confidence levels for replicates

## Key Components

### Java Tools (`src/main/java/org/cruk/ampliconseq/`)

Built with Maven. Four command-line tools using HTSJDK:

- `PileupCounts` - generates pileup summaries from BAM files
- `ExtractAmpliconRegions` - extracts reads for specific amplicon regions
- `AnnotateVcfWithAmpliconIds` - adds amplicon IDs to VCF records
- `AddAssortedAnnotationsToVcf` - adds sequence context and other annotations

Build: `mvn package` (outputs to `target/tools/bin/`)

### R Scripts (`bin/`)

Processing scripts invoked by Nextflow processes:

- `compute_background_noise_thresholds.R` - fits probability distributions for noise filtering
- `apply_background_noise_filters.R` - applies position-level noise filters
- `assess_replicate_vaf.R` - analyzes VAF consistency across replicate libraries
- `summarize_variants.R` - produces final variant table with confidence levels
- `create_qc_report.R` - generates alignment/coverage QC report
- `collate_alignment_coverage_metrics.R` - aggregates Picard metrics
- `create_non_overlapping_amplicon_groups.R` - partitions overlapping amplicons

### Shiny App (`shiny/`)

Interactive visualization for reviewing variant calls. Launch with `start_shiny_app.R`. Provides scatter plots of allele fractions, density plots with fitted distributions, and interactive variant tables.

### Templates (`templates/`)

Shell script templates that define the commands run by each Nextflow process (variant calling, pileup, metrics, annotation, etc.).

## Configuration

Key parameters in `nextflow.config`:

- `variantCaller` - VarDict (default), HaplotypeCaller, or Mutect2
- `minimumAlleleFraction` - minimum allele fraction threshold (default: 0.01)
- `vepAnnotation` - enable Ensembl VEP annotation (default: false)
- `referenceGenomeFasta` - path to indexed reference genome

Execution profiles: `standard` (4 CPU), `bigserver` (40 CPU), `cluster` (Slurm)

## Dependencies

All packaged in Docker container (`docker/`):

- Nextflow 25.10.0+, Java 17+
- GATK 4.6.2.0, VarDict 1.8.3, BCFtools 1.21, Ensembl VEP 115.2
- R 4.3.3 with tidyverse, fitdistrplus, ComplexHeatmap, shiny

## Running

```bash
# Build Java tools
mvn package

# Run pipeline (with container)
nextflow run crukci-bioinformatics/ampliconseq \
    -config ampliconseq.config \
    -with-singularity \
    -profile bigserver
```

Input: sample sheet CSV (ID, Sample, BAM columns) + amplicon coordinates CSV.
Output: `variants.txt`, VCF files, pileup counts, QC report, coverage metrics.
