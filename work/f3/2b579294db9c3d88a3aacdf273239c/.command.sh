#!/bin/bash -euo pipefail
echo "1.6" &> v_pipeline.txt
echo "20.10.0" &> v_nextflow.txt
bismark_genome_preparation --version &> v_bismark_genome_preparation.txt
fastqc --version &> v_fastqc.txt
cutadapt --version &> v_cutadapt.txt
trim_galore --version &> v_trim_galore.txt
bismark --version &> v_bismark.txt
deduplicate_bismark --version &> v_deduplicate_bismark.txt
bismark_methylation_extractor --version &> v_bismark_methylation_extractor.txt
bismark2report --version &> v_bismark2report.txt
bismark2summary --version &> v_bismark2summary.txt
samtools --version &> v_samtools.txt
hisat2 --version &> v_hisat2.txt
bwa &> v_bwa.txt 2>&1 || true
bwameth.py --version &> v_bwameth.txt
picard MarkDuplicates --version &> v_picard_markdups.txt 2>&1 || true
MethylDackel --version &> v_methyldackel.txt
qualimap --version &> v_qualimap.txt || true
preseq &> v_preseq.txt
multiqc --version &> v_multiqc.txt
scrape_software_versions.py &> software_versions_mqc.yaml
