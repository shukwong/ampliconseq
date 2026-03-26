#!/bin/bash

set -e -o pipefail

# merge all control BAM pileup VCFs
bcftools merge \
    --output-type z \
    --output !{merged_pileup_vcf} \
    !{pileup_vcfs}

tabix !{merged_pileup_vcf}

# compute aggregate PON_RefDepth and PON_AltDepth across all control samples
bcftools +fill-tags -Ov !{merged_pileup_vcf} -- -t "PON_RefDepth=sum(RD)" | \
bcftools +fill-tags -Oz -o !{pon_total_counts_vcf} -- -t "PON_AltDepth=sum(AD)"

tabix !{pon_total_counts_vcf}
