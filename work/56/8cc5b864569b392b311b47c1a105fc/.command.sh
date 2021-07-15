#!/bin/bash -euo pipefail
dv_postprocess_variants.py   --ref chr20.fa.gz   --infile call_variants_output.tfrecord   --outfile "NA12878_S1.chr20.10_10p1mb.bam.vcf"
