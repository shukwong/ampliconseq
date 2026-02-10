#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Makes minor modifications to the variant table output by GATK VariantsToTable
# including:
# - add ID column
# - separate allele depths into separate columns
# - rounds quality scores
# - computes the allele fraction from the allelic depths if no values extracted

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--input"), dest = "input_file",
              help = "Variant table from GATK VariantsToTable (AMPLICON, CHROM, POS, REF, ALT, MULTI-ALLELIC, FILTER, QUAL, {ID}.DP {ID}.AD {ID}.AF columns required)"),

  make_option(c("--id"), dest = "id", help = "Library ID"),

  make_option(c("--sample-name"), dest = "sample_name",
              help = "Sample name used in VCF (defaults to ID if not specified)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output file containing all variants including rows for missing calls")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

id <- opt$id
sample_name <- opt$sample_name
input_file <- opt$input_file
output_file <- opt$output_file

if (is.null(id)) stop("Library ID must be specified")
if (is.null(input_file)) stop("Input variants file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

# use sample_name for matching VCF column names; fall back to id if not specified
if (is.null(sample_name)) sample_name <- id

suppressPackageStartupMessages(library(tidyverse))

# read variants file
variants <- read_tsv(input_file, col_types = cols(.default = "c"))

# rename columns using name-based matching rather than positional assignment
# (positional assignment fails if VariantsToTable produces extra columns, e.g.
# when it splits multi-valued genotype fields like AD into separate columns)
dp_col <- str_c(sample_name, ".DP")
ad_col <- str_c(sample_name, ".AD")
af_col <- str_c(sample_name, ".AF")

# check for split AD columns (e.g. {sample}.AD.REF and {sample}.AD.ALT instead of {sample}.AD)
ad_ref_col <- str_c(sample_name, ".AD.REF")
ad_alt_col <- str_c(sample_name, ".AD.ALT")
split_ad <- ad_ref_col %in% colnames(variants) && ad_alt_col %in% colnames(variants)

if (split_ad) {
  variants <- variants %>%
    rename(
      Amplicon = AMPLICON,
      Chromosome = CHROM,
      Position = POS,
      Ref = REF,
      Alt = ALT,
      Multiallelic = MULTIALLELIC,
      Filters = FILTER,
      Quality = QUAL,
      Depth = !!dp_col,
      `Ref depth` = !!ad_ref_col,
      `Alt depth` = !!ad_alt_col,
      `Allele fraction` = !!af_col
    )
} else {
  variants <- variants %>%
    rename(
      Amplicon = AMPLICON,
      Chromosome = CHROM,
      Position = POS,
      Ref = REF,
      Alt = ALT,
      Multiallelic = MULTIALLELIC,
      Filters = FILTER,
      Quality = QUAL,
      Depth = !!dp_col,
      AD = !!ad_col,
      `Allele fraction` = !!af_col
    )
}

# drop any extra columns not needed
expected_cols <- c("Amplicon", "Chromosome", "Position", "Ref", "Alt", "Multiallelic",
                   "Filters", "Quality", "FivePrimeContext", "Depth",
                   if (split_ad) c("Ref depth", "Alt depth") else "AD",
                   "Allele fraction")
variants <- variants %>% select(all_of(expected_cols))

# convert Multiallelic column to boolean and change NAs to FALSE
variants <- variants %>%
  mutate(Multiallelic = as.logical(Multiallelic)) %>%
  mutate(Multiallelic = replace_na(Multiallelic, FALSE))

# separate Ref and Alt allelic depths from AD column (if not already split)
if (!split_ad) {
  variants <- separate(variants, AD, into = c("Ref depth", "Alt depth"), sep = ",")
}

# round quality scores
variants <- mutate(variants, Quality = round(parse_double(Quality)))

# compute allele fraction if all values are missing, i.e. the variant caller
# doesn't create an AF field
if (all(is.na(variants$`Allele fraction`))) {
  variants <- variants %>%
    mutate(across(c(`Ref depth`, `Alt depth`, Depth), parse_integer)) %>%
    mutate(`Allele fraction` = round(`Alt depth` / Depth, digits = 5))
}

# add ID column
variants <- variants %>%
  mutate(ID = id) %>%
  select(ID, everything())

# Fix deletions where the alternate allele is '*' that come from multi-allelic
# variants by adding the preceding base to the reference allele and replacing
# the '*' with the preceding base.
# e.g. chr17 29563076 T * becomes chr17 29563076 GT T
# This requires the FivePrimeContext field to have been set, e.g. by running the
# add-assorted-annotations-to-vcf tool.
variants <- variants %>%
  mutate(Position = as.integer(Position)) %>%
  mutate(Position = ifelse(Alt == "*", Position - 1, Position)) %>%
  mutate(Ref = ifelse(Alt == "*", str_c(FivePrimeContext, Ref), Ref)) %>%
  mutate(Alt = ifelse(Alt == "*", FivePrimeContext, Alt)) %>%
  select(-FivePrimeContext)

# remove duplicates due to splitting multi-allelic variants where they are also
# called separately (have seen this with GATK HaplotypeCaller)
duplicates <- variants %>%
  count(ID, Amplicon, Chromosome, Position, Ref, Alt) %>%
  filter(n > 1)

multiallelic_duplicates <- variants %>%
  semi_join(duplicates, by = c("ID", "Amplicon", "Chromosome", "Ref", "Alt")) %>%
  filter(!Multiallelic)

variants <- variants %>%
  anti_join(multiallelic_duplicates, by = c("ID", "Amplicon", "Chromosome", "Ref", "Alt", "Multiallelic"))

duplicates <- variants %>%
  count(ID, Amplicon, Chromosome, Position, Ref, Alt) %>%
  filter(n > 1)

if (nrow(duplicates) > 0) {
  stop("Duplicates remain in variant table after removing multi-allelic duplicates")
}

# write to output file
write_tsv(variants, output_file, na = "")

