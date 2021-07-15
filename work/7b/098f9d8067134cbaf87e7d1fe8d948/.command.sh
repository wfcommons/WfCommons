#!/bin/bash -euo pipefail
samtools sort SRR389222_sub1_trimmed_bismark_bt2.bam \
    -@ 1  \
    -o SRR389222_sub1_trimmed_bismark_bt2.sorted.bam
preseq lc_extrap -v -B SRR389222_sub1_trimmed_bismark_bt2.sorted.bam -o SRR389222_sub1_trimmed_bismark_bt2.ccurve.txt
