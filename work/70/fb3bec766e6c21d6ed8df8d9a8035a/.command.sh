#!/bin/bash -euo pipefail
gatk --java-options -Xmx6g         ApplyBQSR         -R human_g1k_v37_decoy.small.fasta         --input 9876T.md.bam         --output 3_1-200000_9876T.recal.bam         -L 3_1-200000.bed         --bqsr-recal-file 9876T.recal.table
