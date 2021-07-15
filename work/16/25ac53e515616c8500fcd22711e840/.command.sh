#!/bin/bash -euo pipefail
mkdir -p trimmed
if [[ "false" = "true" && "false" = "false" ]]; then
	#rename files to list results correctly in MultiQC
	if [ "false" = "true" ]; then
		ln -s "2a_S115_L001_R1_001.fastq.gz 2a_S115_L001_R2_001.fastq.gz" "first-trimming_2a_S115_L001_R1_001.fastq.gz 2a_S115_L001_R2_001.fastq.gz"
	else
		ln -s "2a_S115_L001_R1_001.fastq.gz" "first-trimming_2a_S115_L001_R1_001.fastq.gz"
		ln -s "2a_S115_L001_R2_001.fastq.gz" "first-trimming_2a_S115_L001_R2_001.fastq.gz"
	fi

	cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed 					-o second-trimming_2a_S115_L001_R1_001.fastq.gz -p second-trimming_2a_S115_L001_R2_001.fastq.gz 					first-trimming_2a_S115_L001_R1_001.fastq.gz first-trimming_2a_S115_L001_R2_001.fastq.gz >> cutadapt_log_2a.txt
	cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-trimmed 					-o trimmed/2a_S115_L001_R1_001.fastq.gz -p trimmed/2a_S115_L001_R2_001.fastq.gz 					second-trimming_2a_S115_L001_R1_001.fastq.gz second-trimming_2a_S115_L001_R2_001.fastq.gz >> cutadapt_log_2a.txt
else
	cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed 					-o trimmed/2a_S115_L001_R1_001.fastq.gz -p trimmed/2a_S115_L001_R2_001.fastq.gz 2a_S115_L001_R1_001.fastq.gz 2a_S115_L001_R2_001.fastq.gz 					>> cutadapt_log_2a.txt
fi
