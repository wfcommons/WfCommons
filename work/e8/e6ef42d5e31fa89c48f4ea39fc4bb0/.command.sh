#!/bin/bash -euo pipefail
mkdir ready
[[ `samtools view -H NA12878_S1.chr20.10_10p1mb.bam | grep '@RG' | wc -l`   > 0 ]] && { mv NA12878_S1.chr20.10_10p1mb.bam ready;}|| { picard AddOrReplaceReadGroups   I=NA12878_S1.chr20.10_10p1mb.bam   O=ready/NA12878_S1.chr20.10_10p1mb.bam   RGID=4   RGLB=lib1   RGPL=illumina   RGPU=unit1   RGSM=20;}
cd ready ;samtools index NA12878_S1.chr20.10_10p1mb.bam;
