#!/usr/bin/env Rscript

# Aggregates pileup counts from control BAMs (Panel of Normals) and adds
# PON Depth, PON Alt depth and PON Alt fraction columns to the variants table

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--variants"), dest = "variants_file",
              help = "TSV file containing variants (Amplicon, Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--pon-pileup-counts"), dest = "pon_pileup_counts_file",
              help = "TSV file containing pileup counts from control BAMs. Supports either per-base columns ('Reference base', 'A count'...'T count') or direct allele columns (Ref, Alt, 'Alt count')."),

  make_option(c("--output"), dest = "output_file",
              help = "Output variants file with PON columns added")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

variants_file <- opt$variants_file
pon_pileup_counts_file <- opt$pon_pileup_counts_file
output_file <- opt$output_file

if (is.null(variants_file)) stop("Input variant file must be specified")
if (is.null(pon_pileup_counts_file)) stop("PON pileup counts file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

variants <- read_tsv(variants_file, col_types = cols(.default = "c"))

pon_pileup_counts <- read_tsv(pon_pileup_counts_file, col_types = cols(.default = "c"))

# If Amplicon is missing (variant-targeted pileup), map it from variants
if (!"Amplicon" %in% colnames(pon_pileup_counts)) {
  variants_for_amplicon <- variants %>% distinct(Amplicon, Chromosome, Position, Ref, Alt)
  pon_pileup_counts <- pon_pileup_counts %>%
    left_join(variants_for_amplicon, by = c("Chromosome", "Position", "Ref", "Alt")) %>%
    distinct(ID, Amplicon, Chromosome, Position, Ref, Alt, .keep_all = TRUE) %>%
    relocate(Amplicon, .before = Chromosome)
}

if (all(c("Ref", "Alt", "Alt count") %in% colnames(pon_pileup_counts))) {
  # Variant-targeted pileup already carries the queried allele, so this path
  # works for both SNVs and indels.
  pon_summary <- pon_pileup_counts %>%
    semi_join(variants, by = c("Amplicon", "Chromosome", "Position", "Ref", "Alt")) %>%
    mutate(across(c(Depth, `Alt count`), parse_integer)) %>%
    group_by(Amplicon, Chromosome, Position, Ref, Alt) %>%
    summarize(
      `PON Depth` = sum(Depth, na.rm = TRUE),
      `PON Alt depth` = sum(`Alt count`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(`PON Alt fraction` = `PON Alt depth` / `PON Depth`)
} else {
  # Legacy amplicon-wide pileup path only has nucleotide counts.
  pon_summary <- pon_pileup_counts %>%
    semi_join(variants, by = c("Amplicon", "Chromosome", "Position")) %>%
    select(ID, Amplicon, Chromosome, Position, Ref = `Reference base`, `A count`:`T count`, Depth) %>%
    pivot_longer(`A count`:`T count`, names_to = "Alt", values_to = "Count") %>%
    mutate(Alt = str_remove(Alt, " count$")) %>%
    mutate(across(c(Count, Depth), parse_integer)) %>%
    group_by(Amplicon, Chromosome, Position, Ref, Alt) %>%
    summarize(
      `PON Depth` = sum(Depth, na.rm = TRUE),
      `PON Alt depth` = sum(Count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(`PON Alt fraction` = `PON Alt depth` / `PON Depth`)
}

variants %>%
  left_join(pon_summary, by = c("Amplicon", "Chromosome", "Position", "Ref", "Alt")) %>%
  write_tsv(output_file, na = "")
