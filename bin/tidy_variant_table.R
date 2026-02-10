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

  make_option(c("--output"), dest = "output_file",
              help = "Output file containing all variants including rows for missing calls")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

id <- opt$id
input_file <- opt$input_file
output_file <- opt$output_file

if (is.null(id)) stop("Library ID must be specified")
if (is.null(input_file)) stop("Input variants file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# read variants file
variants <- read_tsv(input_file, col_types = cols(.default = "c"))

# auto-detect sample genotype columns from VariantsToTable output
# the fixed INFO columns are known; genotype columns follow the pattern {sample}.DP, {sample}.AD, {sample}.AF
# when there are multiple samples in the VCF we take the first sample's columns
col_names <- colnames(variants)
dp_cols <- col_names[str_ends(col_names, fixed(".DP"))]
ad_cols <- col_names[str_ends(col_names, fixed(".AD"))]
af_cols <- col_names[str_ends(col_names, fixed(".AF"))]

if (length(dp_cols) == 0) stop("No .DP genotype column found in input")
if (length(af_cols) == 0) stop("No .AF genotype column found in input")

# use the first sample's columns
dp_col <- dp_cols[1]
af_col <- af_cols[1]

# check for split AD columns (e.g. {sample}.AD.REF and {sample}.AD.ALT instead of {sample}.AD)
ad_ref_cols <- col_names[str_ends(col_names, fixed(".AD.REF"))]
ad_alt_cols <- col_names[str_ends(col_names, fixed(".AD.ALT"))]
split_ad <- length(ad_ref_cols) > 0 && length(ad_alt_cols) > 0

if (split_ad) {
  ad_ref_col <- ad_ref_cols[1]
  ad_alt_col <- ad_alt_cols[1]
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
  if (length(ad_cols) == 0) stop("No .AD genotype column found in input")
  ad_col <- ad_cols[1]
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

# drop any extra columns (e.g. from additional samples in multi-sample VCFs)
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
