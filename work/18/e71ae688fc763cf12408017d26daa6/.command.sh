#!/bin/bash -euo pipefail
mkdir -p trimmed
if [[ "false" = "true" && "false" = "false" ]]; then
	#rename files to list results correctly in MultiQC
	if [ "false" = "true" ]; then
		ln -s "1_S103_L001_R1_001.fastq.gz 1_S103_L001_R2_001.fastq.gz" "first-trimming_1_S103_L001_R1_001.fastq.gz 1_S103_L001_R2_001.fastq.gz"
	else
		ln -s "1_S103_L001_R1_001.fastq.gz" "first-trimming_1_S103_L001_R1_001.fastq.gz"
		ln -s "1_S103_L001_R2_001.fastq.gz" "first-trimming_1_S103_L001_R2_001.fastq.gz"
	fi

	cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed 					-o second-trimming_1_S103_L001_R1_001.fastq.gz -p second-trimming_1_S103_L001_R2_001.fastq.gz 					first-trimming_1_S103_L001_R1_001.fastq.gz first-trimming_1_S103_L001_R2_001.fastq.gz >> cutadapt_log_1.txt
	cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-trimmed 					-o trimmed/1_S103_L001_R1_001.fastq.gz -p trimmed/1_S103_L001_R2_001.fastq.gz 					second-trimming_1_S103_L001_R1_001.fastq.gz second-trimming_1_S103_L001_R2_001.fastq.gz >> cutadapt_log_1.txt
else
	cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed 					-o trimmed/1_S103_L001_R1_001.fastq.gz -p trimmed/1_S103_L001_R2_001.fastq.gz 1_S103_L001_R1_001.fastq.gz 1_S103_L001_R2_001.fastq.gz 					>> cutadapt_log_1.txt
fi
