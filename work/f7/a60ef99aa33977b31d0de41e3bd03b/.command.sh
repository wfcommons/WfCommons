#!/bin/bash -euo pipefail
fastqc -t 2 -q 1234N_1234N_M6_R1.fastq.gz 1234N_1234N_M6_R2.fastq.gz
