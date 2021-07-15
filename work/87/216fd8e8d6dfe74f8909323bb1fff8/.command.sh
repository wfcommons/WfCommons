#!/bin/bash -euo pipefail
samtools sort SRR389222_sub2_trimmed_bismark_bt2.deduplicated.bam \
    -@ 2  \
    -o SRR389222_sub2_trimmed_bismark_bt2.deduplicated.sorted.bam
qualimap bamqc  \
    -bam SRR389222_sub2_trimmed_bismark_bt2.deduplicated.sorted.bam \
    -outdir SRR389222_sub2_trimmed_bismark_bt2.deduplicated_qualimap \
    --collect-overlap-pairs \
    --java-mem-size=6G \
    -nt 2
