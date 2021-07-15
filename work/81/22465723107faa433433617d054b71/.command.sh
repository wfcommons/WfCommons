#!/bin/bash -euo pipefail
multiqc -f    . \
    -m custom_content -m picard -m qualimap -m bismark -m samtools -m preseq -m cutadapt -m fastqc
