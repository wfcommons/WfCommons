#!/bin/bash -euo pipefail
qualimap --java-mem-size=6G         bamqc         -bam 1234N_1234N_M1.bam         --paint-chromosome-limits         --genome-gc-distr HUMAN                  -nt 2         -skip-duplicated         --skip-dup-mode 0         -outdir 1234N_1234N_M1         -outformat HTML
