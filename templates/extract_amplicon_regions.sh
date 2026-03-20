#!/bin/bash

set -e -o pipefail

EXTRA_ARGS=""
if [ !{params.maximumReadsPerAmplicon} -gt 0 ]; then
    EXTRA_ARGS="${EXTRA_ARGS} --maximum-reads-per-amplicon !{params.maximumReadsPerAmplicon}"
fi
if [ !{params.downsamplePercentile} -gt 0 ]; then
    EXTRA_ARGS="${EXTRA_ARGS} --downsample-percentile !{params.downsamplePercentile}"
fi

JAVA_OPTS="-Xmx!{java_mem}m" extract-amplicon-regions \
    --id "!{id}" \
    --input !{bam} \
    --amplicon-intervals !{amplicon_bed} \
    --output "!{amplicon_bam}" \
    --coverage "!{amplicon_coverage}" \
    --maximum-distance !{params.maxDistanceFromAmpliconEnd} \
    --require-both-ends-anchored=!{params.requireBothEndsAnchored} \
    --unmark-duplicate-reads \
    ${EXTRA_ARGS}
