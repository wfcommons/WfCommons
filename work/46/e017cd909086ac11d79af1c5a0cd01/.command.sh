#!/bin/bash -euo pipefail
fastqc -t 2 -q 1234N_1234N_M4_R1.fastq.gz 1234N_1234N_M4_R2.fastq.gz
