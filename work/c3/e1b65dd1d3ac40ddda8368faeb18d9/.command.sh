#!/bin/bash -euo pipefail
dv_call_variants.py     --cores 2     --sample NA12878_S1.chr20.10_10p1mb.bam     --outfile NA12878_S1.chr20.10_10p1mb_call_variants_output.tfrecord     --examples NA12878_S1.chr20.10_10p1mb_shardedExamples     --model wgs
