#!/bin/bash -euo pipefail
fastqc -t 2 -q 9876T_9876T_M4_R1.fastq.gz 9876T_9876T_M4_R2.fastq.gz
