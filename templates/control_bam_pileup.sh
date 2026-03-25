#!/bin/bash

set -e -o pipefail

# use sample name provided by the pipeline
sample_name="!{sample_name}"

# check if the fake VCF has any variants
variant_count=$(zgrep -v '^#' !{fake_vcf} | wc -l || true)

if [ "${variant_count}" -lt 1 ]; then
    # create an empty pileup VCF with proper headers
    printf "##fileformat=VCFv4.2\n" > pileup.vcf
    printf "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n" >> pileup.vcf
    printf "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth matching reference\">\n" >> pileup.vcf
    printf "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth matching alternate\">\n" >> pileup.vcf
    printf "##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant frequency\">\n" >> pileup.vcf
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_name}\n" >> pileup.vcf
else
    # decompress fake VCF if needed for GetBaseCountsMultiSample
    fake_vcf_unzipped=fake_input.vcf
    zcat !{fake_vcf} > ${fake_vcf_unzipped}

    /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample \
        --fasta !{reference_sequence} \
        --bam ${sample_name}:!{bam} \
        --vcf ${fake_vcf_unzipped} \
        --output pileup.vcf \
        --maq !{params.ponPileupMappingQuality} \
        --baq !{params.ponPileupBaseQuality} \
        --thread !{task.cpus} \
        --max_block_dist 10000
fi

bgzip pileup.vcf && tabix pileup.vcf.gz
mv pileup.vcf.gz !{pileup_vcf}
mv pileup.vcf.gz.tbi !{pileup_vcf_tbi}
