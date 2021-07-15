#!/bin/bash -euo pipefail
gatk --java-options "-Xms3g -Xmx5g"         MarkDuplicates         --INPUT 9876T.bam         --METRICS_FILE 9876T.bam.metrics         --TMP_DIR .         --ASSUME_SORT_ORDER coordinate         --CREATE_INDEX true         --OUTPUT 9876T.md.bam

mv 9876T.md.bai 9876T.md.bam.bai
