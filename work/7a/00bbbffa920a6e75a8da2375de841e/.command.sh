#!/bin/bash -euo pipefail
gatk --java-options -Xmx6g         BaseRecalibrator         -I 9876T.md.bam         -O 8_1-1276_9876T.recal.table         --tmp-dir .         -R human_g1k_v37_decoy.small.fasta         -L 8_1-1276.bed         --known-sites dbsnp_138.b37.small.vcf.gz         --known-sites Mills_1000G_gold_standard_and_1000G_phase1.indels.b37.small.vcf.gz         --verbosity INFO
