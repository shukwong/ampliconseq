#!/bin/bash

set -e -o pipefail

# Get the actual file size in bytes, following symlinks
actual_file=$(readlink -f !{bam})
file_size=$(stat -L -c %s "${actual_file}")

# Convert max size from GB to bytes using awk
max_size_bytes=$(awk "BEGIN {printf \"%.0f\", !{max_size_gb} * 1024 * 1024 * 1024}")

if [ "${file_size}" -gt "${max_size_bytes}" ]; then
    echo "BAM file size (${file_size} bytes) exceeds maximum (!{max_size_gb}GB). Subsampling..."

    # Calculate the sampling fraction to achieve target size using awk
    fraction=$(awk "BEGIN {printf \"%.6f\", ${max_size_bytes} / ${file_size}}")

    # Use a random seed based on the sample ID for reproducibility
    seed=$(echo "!{id}" | cksum | cut -d' ' -f1)
    seed=$((seed % 10000))

    # Subsample using samtools view with the calculated fraction
    # Format: seed.fraction (e.g., 42.123456)
    samtools view -b -s ${seed}${fraction#0} -o !{output_bam} !{bam}

    # Index the subsampled BAM
    samtools index !{output_bam}

    echo "Subsampling complete. New file size: $(stat -c %s !{output_bam}) bytes"
else
    echo "BAM file size (${file_size} bytes) is within limit (!{max_size_gb}GB). Creating symlink..."

    # Just create a symlink to the original file
    ln -s "${actual_file}" !{output_bam}

    # Also link the index if it exists, otherwise create one
    if [ -f "${actual_file}.bai" ]; then
        ln -s "${actual_file}.bai" !{output_bam}.bai
    elif [ -f "${actual_file%.bam}.bai" ]; then
        ln -s "${actual_file%.bam}.bai" !{output_bam}.bai
    else
        samtools index !{output_bam}
    fi
fi
