#!/bin/bash -euo pipefail
bismark SRR389222_sub2_trimmed.fq.gz \
    --bowtie2 \
    --bam        \
    --genome BismarkIndex \
    SRR389222_sub2_trimmed.fq.gz \
     \
