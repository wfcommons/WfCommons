#!/bin/bash -euo pipefail
samtools merge --threads 2 1234N.recal.bam 3_1-200000_1234N.recal.bam 11_1-3679_1234N.recal.bam 1_1-200000_1234N.recal.bam 8_1-1276_1234N.recal.bam 2_1-200000_1234N.recal.bam X_1-200000_1234N.recal.bam
samtools index 1234N.recal.bam
