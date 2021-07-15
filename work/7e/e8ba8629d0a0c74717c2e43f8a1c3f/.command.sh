#!/bin/bash -euo pipefail
echo 1.0 &> v_nf_deepvariant.txt
echo 20.10.0 &> v_nextflow.txt
ls /opt/conda/pkgs/ &> v_deepvariant.txt
python --version &> v_python.txt
pip --version &> v_pip.txt
samtools --version &> v_samtools.txt
lbzip2 --version &> v_lbzip2.txt
bzip2 --version &> v_bzip2.txt
scrape_software_versions.py &> software_versions_mqc.yaml
