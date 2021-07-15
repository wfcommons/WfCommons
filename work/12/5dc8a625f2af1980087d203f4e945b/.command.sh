#!/bin/bash -euo pipefail
gatk --java-options -Xmx6g         GatherBQSRReports         -I 11_1-3679_1234N.recal.table -I 3_1-200000_1234N.recal.table -I 1_1-200000_1234N.recal.table -I X_1-200000_1234N.recal.table -I 2_1-200000_1234N.recal.table -I 8_1-1276_1234N.recal.table         -O 1234N.recal.table
