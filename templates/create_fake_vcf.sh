#!/bin/bash

set -e -o pipefail

# create a simplified "fake" VCF containing only variant positions
# suitable for GetBaseCountsMultiSample pileup
echo -e "##fileformat=VCFv4.2" > fake.vcf
echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> fake.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> fake.vcf

zcat !{merged_vcf} \
    | grep -v '^#' \
    | awk '{print $1, $2, $3, $4, $5, $6, "PASS\t.\tGT\t0/1"}' OFS='\t' \
    >> fake.vcf

bgzip fake.vcf && tabix fake.vcf.gz
mv fake.vcf.gz !{fake_vcf}
mv fake.vcf.gz.tbi !{fake_vcf_tbi}
