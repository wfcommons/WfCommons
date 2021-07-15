#!/bin/bash -euo pipefail
gatk --java-options "-Xmx6g"         CreateSequenceDictionary         --REFERENCE human_g1k_v37_decoy.small.fasta         --OUTPUT human_g1k_v37_decoy.small.dict
