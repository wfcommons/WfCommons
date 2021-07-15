#!/bin/bash -euo pipefail
samtools merge --threads 2 9876T.recal.bam 1_1-200000_9876T.recal.bam 11_1-3679_9876T.recal.bam 2_1-200000_9876T.recal.bam 8_1-1276_9876T.recal.bam 3_1-200000_9876T.recal.bam X_1-200000_9876T.recal.bam
samtools index 9876T.recal.bam
