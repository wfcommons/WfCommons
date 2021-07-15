#!/bin/bash -euo pipefail
echo 1.2.2 > v_pipeline.txt
echo 20.10.0 > v_nextflow.txt
bowtie2 --version > v_bowtie2.txt
python --version > v_python.txt 2>&1
samtools --version > v_samtools.txt
multiqc --version > v_multiqc.txt
scrape_software_versions.py &> software_versions_mqc.yaml
