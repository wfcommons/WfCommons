#!/bin/bash -euo pipefail
qualimap --java-mem-size=6G         bamqc         -bam 9876T_9876T_M3.bam         --paint-chromosome-limits         --genome-gc-distr HUMAN                  -nt 2         -skip-duplicated         --skip-dup-mode 0         -outdir 9876T_9876T_M3         -outformat HTML
