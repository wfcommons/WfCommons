#!/bin/bash -euo pipefail
bismark2report \
    --alignment_report SRR389222_sub2_trimmed_bismark_bt2_SE_report.txt \
    --dedup_report SRR389222_sub2_trimmed_bismark_bt2.deduplication_report.txt \
    --splitting_report SRR389222_sub2_trimmed_bismark_bt2.deduplicated_splitting_report.txt \
    --mbias_report SRR389222_sub2_trimmed_bismark_bt2.deduplicated.M-bias.txt
