#!/bin/bash -euo pipefail
bismark_methylation_extractor   \
       \
    --bedGraph \
    --counts \
    --gzip \
    -s \
    --report \
    SRR389222_sub3_trimmed_bismark_bt2.deduplicated.bam
