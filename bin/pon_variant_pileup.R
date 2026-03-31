#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------
# Extract per-allele counts from bcftools pileup/call output and map Amplicon
# from the variants table so we can reuse add_pon_pileup.R for aggregation.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--variants"), dest = "variants_file",
              help = "TSV/CSV file containing variants (Amplicon, Chromosome, Position, Ref, Alt required)"),
  make_option(c("--pileup"), dest = "pileup_file",
              help = "Raw pileup TSV from bcftools query with columns: CHROM POS REF ALT SAMPLE DP AD"),
  make_option(c("--output"), dest = "output_file",
              help = "Output TSV with columns: ID, Amplicon, Chromosome, Position, Ref, Alt, Depth, Alt count")
)

opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

variants_file <- opt$variants_file
pileup_file <- opt$pileup_file
output_file <- opt$output_file

if (is.null(variants_file)) stop("Variants file must be specified")
if (is.null(pileup_file)) stop("Pileup file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# load variants (needed only for Amplicon mapping)
variants <- if (str_ends(str_to_lower(variants_file), "\\.csv")) {
  read_csv(variants_file, col_types = cols(.default = col_character()))
} else {
  read_tsv(variants_file, col_types = cols(.default = col_character()))
}

variants_subset <- variants %>%
  distinct(Amplicon, Chromosome, Position, Ref, Alt)

# load pileup query output
col_names <- c("Chromosome", "Position", "Ref", "Alt", "ID", "Depth", "AD")
pileup <- read_tsv(pileup_file, col_names = col_names, col_types = cols(.default = "c"))

# split AD (ref,alt,alt2...) and keep only the current ALT count; the VCFs
# were normalized before querying so each row should be biallelic.
split_counts <- function(ad_string, depth_string) {
  ads <- str_split(ad_string, ",", simplify = TRUE)
  depth <- suppressWarnings(as.integer(depth_string))
  alt_count <- suppressWarnings(as.integer(ads[2]))
  list(depth = depth, alt_count = alt_count)
}

pileup_parsed <- pileup %>%
  rowwise() %>%
  mutate(parsed = list(split_counts(AD, Depth))) %>%
  mutate(Depth = parsed$depth,
         `Alt count` = parsed$alt_count,
         .keep = "unused") %>%
  ungroup() %>%
  select(-parsed)

# add Amplicon from variants table
pileup_with_amplicon <- pileup_parsed %>%
  left_join(variants_subset, by = c("Chromosome", "Position", "Ref", "Alt")) %>%
  select(ID, Amplicon, Chromosome, Position, Ref, Alt, Depth, `Alt count`)

write_tsv(pileup_with_amplicon, output_file, na = "")
