#!/bin/bash -euo pipefail
echo 1.2.0 > v_pipeline.txt
echo 20.10.0 > v_nextflow.txt
fastqc --version > v_fastqc.txt
multiqc --version > v_multiqc.txt
cutadapt --version > v_cutadapt.txt
qiime --version > v_qiime.txt
scrape_software_versions.py &> software_versions_mqc.yaml
