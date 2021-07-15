#!/bin/bash -euo pipefail
gatk --java-options -Xmx6g         GatherBQSRReports         -I X_1-200000_9876T.recal.table -I 2_1-200000_9876T.recal.table -I 8_1-1276_9876T.recal.table -I 3_1-200000_9876T.recal.table -I 1_1-200000_9876T.recal.table -I 11_1-3679_9876T.recal.table         -O 9876T.recal.table
