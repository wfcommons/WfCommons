#!/bin/bash -euo pipefail
samtools view -@ 2 -h -f 0x40 example_pe.bam > output_1.bam
samtools view -@ 2 -h -f 0x80 example_pe.bam > output_2.bam
samtools bam2fq output_1.bam > output_1.fastq
samtools bam2fq output_2.bam > output_2.fastq
yara_mapper -e 3 -t 2 -f bam yara/hla_reference_dna output_1.fastq output_2.fastq > output.bam
samtools view -@ 2 -h -F 4 -f 0x40 -b1 output.bam > mapped_1.bam
samtools view -@ 2 -h -F 4 -f 0x80 -b1 output.bam > mapped_2.bam
