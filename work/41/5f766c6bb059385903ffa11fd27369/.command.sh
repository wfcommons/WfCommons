#!/bin/bash -euo pipefail
alleleCounter --version &> v_allelecount.txt 2>&1 || true
bcftools --version &> v_bcftools.txt 2>&1 || true
bwa version &> v_bwa.txt 2>&1 || true
cnvkit.py version &> v_cnvkit.txt 2>&1 || true
configManta.py --version &> v_manta.txt 2>&1 || true
configureStrelkaGermlineWorkflow.py --version &> v_strelka.txt 2>&1 || true
echo "2.7" &> v_pipeline.txt 2>&1 || true
echo "20.10.0" &> v_nextflow.txt 2>&1 || true
snpEff -version &> v_snpeff.txt 2>&1 || true
fastqc --version &> v_fastqc.txt 2>&1 || true
freebayes --version &> v_freebayes.txt 2>&1 || true
freec &> v_controlfreec.txt 2>&1 || true
gatk ApplyBQSR --help &> v_gatk.txt 2>&1 || true
msisensor &> v_msisensor.txt 2>&1 || true
multiqc --version &> v_multiqc.txt 2>&1 || true
qualimap --version &> v_qualimap.txt 2>&1 || true
R --version &> v_r.txt 2>&1 || true
R -e "library(ASCAT); help(package='ASCAT')" &> v_ascat.txt 2>&1 || true
samtools --version &> v_samtools.txt 2>&1 || true
tiddit &> v_tiddit.txt 2>&1 || true
trim_galore -v &> v_trim_galore.txt 2>&1 || true
vcftools --version &> v_vcftools.txt 2>&1 || true
vep --help &> v_vep.txt 2>&1 || true

scrape_software_versions.py &> software_versions_mqc.yaml
