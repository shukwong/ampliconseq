#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------
# Extract per-allele counts from raw bcftools mpileup output and produce a
# per-(control-sample, variant) table that add_pon_pileup.R can aggregate.
#
# Raw mpileup discovers ALT alleles from reads, so a queried variant's ALT may
# not appear if the control BAM has 0 supporting reads.  We handle this by
# left-joining from the variants table to the pileup on (Chromosome, Position)
# for depth, and on (Chromosome, Position, Alt) for alt-specific counts,
# defaulting to 0 when no matching ALT is found.
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

# load variants
variants <- if (str_ends(str_to_lower(variants_file), "\\.csv")) {
  read_csv(variants_file, col_types = cols(.default = col_character()))
} else {
  read_tsv(variants_file, col_types = cols(.default = col_character()))
}

variants_subset <- variants %>%
  distinct(Amplicon, Chromosome, Position, Ref, Alt)

# load pileup query output from bcftools query
col_names <- c("Chromosome", "Position", "Ref", "Alt", "ID", "Depth", "AD")
pileup <- read_tsv(pileup_file, col_names = col_names, col_types = cols(.default = "c"))

# -- per-position depth table (one row per sample per position) ----------------
pileup_depth <- pileup %>%
  mutate(Depth = suppressWarnings(as.integer(Depth))) %>%
  distinct(ID, Chromosome, Position, .keep_all = TRUE) %>%
  select(ID, Chromosome, Position, Depth)

# -- per-allele alt count table ------------------------------------------------
# Expand multiallelic rows (ALT="G,T", AD="100,3,2") into one row per ALT allele
# with the corresponding alt count from the AD field.
expand_alleles <- function(alt_string, ad_string) {
  alts <- str_split(alt_string, ",")[[1]]
  ads  <- str_split(ad_string, ",")[[1]]
  # AD format: ref_count, alt1_count, alt2_count, ...
  # Skip ads[1] (ref count); pair remaining with ALTs
  alt_counts <- suppressWarnings(as.integer(ads[-1]))
  n <- min(length(alts), length(alt_counts))
  if (n == 0 || all(alts == ".")) {
    return(tibble(Alt = character(0), `Alt count` = integer(0)))
  }
  tibble(Alt = alts[1:n], `Alt count` = alt_counts[1:n])
}

pileup_alleles <- pileup %>%
  filter(!is.na(Alt), Alt != ".") %>%
  rowwise() %>%
  mutate(expanded = list(expand_alleles(Alt, AD))) %>%
  ungroup() %>%
  unnest(expanded) %>%
  select(ID, Chromosome, Position, Alt, `Alt count`)

# -- build output: one row per (control sample, queried variant) ---------------
control_samples <- pileup %>% distinct(ID) %>% pull(ID)

result <- variants_subset %>%
  crossing(ID = control_samples) %>%
  # join position-level depth
  left_join(pileup_depth, by = c("ID", "Chromosome", "Position")) %>%
  # join allele-specific count
  left_join(pileup_alleles, by = c("ID", "Chromosome", "Position", "Alt")) %>%
  # default missing alt count to 0; missing depth to 0
  mutate(
    `Alt count` = replace_na(`Alt count`, 0L),
    Depth = replace_na(Depth, 0L)
  ) %>%
  select(ID, Amplicon, Chromosome, Position, Ref, Alt, Depth, `Alt count`)

write_tsv(result, output_file, na = "")
