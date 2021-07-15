#!/bin/bash -euo pipefail
gatk --java-options -Xmx6g         ApplyBQSR         -R human_g1k_v37_decoy.small.fasta         --input 1234N.md.bam         --output X_1-200000_1234N.recal.bam         -L X_1-200000.bed         --bqsr-recal-file 1234N.recal.table
