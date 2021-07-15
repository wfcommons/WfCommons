#!/bin/bash -euo pipefail
echo 1.2.0 > v_pipeline.txt
echo 20.10.0 > v_nextflow.txt
multiqc --version &> v_multiqc.txt 2>&1 || true
samtools --version &> v_samtools.txt 2>&1 || true
yara_mapper --help  &> v_yara.txt 2>&1 || true
cat $(which OptiTypePipeline.py) &> v_optitype.txt 2>&1 || true
scrape_software_versions.py &> software_versions_mqc.yaml
