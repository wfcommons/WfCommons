#!/bin/bash -euo pipefail
mkdir logs
mkdir NA12878_S1.chr20.10_10p1mb_shardedExamples
dv_make_examples.py   --cores 2   --sample NA12878_S1.chr20.10_10p1mb.bam   --ref chr20.fa.gz   --reads NA12878_S1.chr20.10_10p1mb.bam   --regions test_nist.b37_chr20_100kbp_at_10mb.bed   --logdir logs   --examples NA12878_S1.chr20.10_10p1mb_shardedExamples
