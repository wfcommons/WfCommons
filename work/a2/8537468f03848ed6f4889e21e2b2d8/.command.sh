#!/bin/bash -euo pipefail
gatk --java-options -Xmx6g         ApplyBQSR         -R human_g1k_v37_decoy.small.fasta         --input 9876T.md.bam         --output 8_1-1276_9876T.recal.bam         -L 8_1-1276.bed         --bqsr-recal-file 9876T.recal.table
