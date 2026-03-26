#!/bin/bash

set -e -o pipefail

# bgzip and index each VCF, then merge all per-library VCFs into one
for vcf in !{vcfs}
do
    bgzip -c ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    echo ${vcf}.gz >> vcf_list.txt
done

# merge all sample VCFs into a single multi-sample VCF
bcftools merge \
    --file-list vcf_list.txt \
    --output-type z \
    --output !{merged_vcf}

tabix !{merged_vcf}
