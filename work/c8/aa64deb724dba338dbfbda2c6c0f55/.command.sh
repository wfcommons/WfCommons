#!/bin/bash -euo pipefail
echo 1.0.2 > v_pipeline.txt
echo 20.10.0 > v_nextflow.txt
fastqc --version > v_fastqc.txt
multiqc --version > v_multiqc.txt
STAR --version > v_star.txt
bowtie --version > v_bowtie.txt
cutadapt --version > v_cutadapt.txt
samtools --version > v_samtools.txt
bedtools --version > v_bedtools.txt
read_distribution.py --version > v_rseqc.txt
sortmerna --version > v_sortmerna.txt
scrape_software_versions.py &> software_versions_mqc.yaml
