#!/bin/bash
for i in {1..9}; do echo "======= $i"; time /home/tgcoleman/1000genome-sequential/bin/individuals.py ALL.chr1.250000.vcf 1 1 1001 3000 -s /home/tgcoleman/wfcommons/wfcommons/wfperf/new_test/real_individuals_$i; done
