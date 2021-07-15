#!/bin/bash -euo pipefail
bwa mem -K 100000000 -R "@RG\tID:9876T_M3\tPU:9876T_M3\tSM:9876T\tLB:9876T\tPL:illumina" -B 3 -t 2 -M human_g1k_v37_decoy.small.fasta     C09D_T_111207.6.ACCAACTG_R1_xxx.fastq.gz C09D_T_111207.6.ACCAACTG_R2_xxx.fastq.gz |     samtools sort --threads 2 -m 2G - > 9876T_9876T_M3.bam
