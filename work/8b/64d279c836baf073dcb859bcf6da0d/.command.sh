#!/bin/bash -euo pipefail
gatk --java-options "-Xms3g -Xmx5g"         MarkDuplicates         --INPUT 1234N.bam         --METRICS_FILE 1234N.bam.metrics         --TMP_DIR .         --ASSUME_SORT_ORDER coordinate         --CREATE_INDEX true         --OUTPUT 1234N.md.bam

mv 1234N.md.bai 1234N.md.bam.bai
