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
              help = "Output TSV with columns: ID, Amplicon, Chromosome, Position, Reference base, Depth, A count, C count, G count, T count")
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

# split AD (ref,alt,alt2...) and compute counts per nucleotide
split_counts <- function(ref_base, alt_base, ad_string, depth_string) {
  ads <- str_split(ad_string, ",", simplify = TRUE)
  depth <- suppressWarnings(as.integer(depth_string))
  ref_count <- suppressWarnings(as.integer(ads[1]))
  alt_count <- suppressWarnings(as.integer(ads[2]))
  counts <- c(A = 0L, C = 0L, G = 0L, T = 0L)
  if (!is.na(ref_count) && ref_base %in% names(counts)) counts[ref_base] <- counts[ref_base] + ref_count
  if (!is.na(alt_count) && alt_base %in% names(counts)) counts[alt_base] <- counts[alt_base] + alt_count
  list(depth = depth, counts = counts)
}

pileup_parsed <- pileup %>%
  rowwise() %>%
  mutate(parsed = list(split_counts(Ref, Alt, AD, Depth))) %>%
  mutate(Depth = parsed$depth,
         `A count` = parsed$counts["A"],
         `C count` = parsed$counts["C"],
         `G count` = parsed$counts["G"],
         `T count` = parsed$counts["T"],
         .keep = "unused") %>%
  ungroup() %>%
  select(-parsed)

# add Amplicon from variants table
pileup_with_amplicon <- pileup_parsed %>%
  left_join(variants_subset, by = c("Chromosome", "Position", "Ref", "Alt")) %>%
  rename(`Reference base` = Ref)

write_tsv(pileup_with_amplicon, output_file, na = "")

