#!/bin/bash -euo pipefail
bwa mem -K 100000000 -R "@RG\tID:1234N_M2\tPU:1234N_M2\tSM:1234N\tLB:1234N\tPL:illumina"  -t 2 -M human_g1k_v37_decoy.small.fasta     C097F_N_111207.2.AGTTGCTT_R1_xxx.fastq.gz C097F_N_111207.2.AGTTGCTT_R2_xxx.fastq.gz |     samtools sort --threads 2 -m 2G - > 1234N_1234N_M2.bam
