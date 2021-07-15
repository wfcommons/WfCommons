#!/bin/bash -euo pipefail
export HOME="${PWD}/HOME"
IFS=',' read -r -a trunclen <<< "230,229"

#denoise samples with DADA2 and produce
qiime dada2 denoise-paired  			--i-demultiplexed-seqs demux.qza  			--p-trunc-len-f ${trunclen[0]} 			--p-trunc-len-r ${trunclen[1]} 			--p-max-ee-f 2 			--p-max-ee-r 2 			--p-n-threads 0  			--o-table table.qza  			--o-representative-sequences rep-seqs.qza  			--o-denoising-stats stats.qza 			--verbose 		>dada_report.txt

#produce dada2 stats "dada_stats/stats.tsv"
qiime tools export --input-path stats.qza 			--output-path dada_stats

#produce raw count table in biom format "table/feature-table.biom"
qiime tools export --input-path table.qza  			--output-path table

#produce raw count table
biom convert -i table/feature-table.biom 			-o table/feature-table.tsv  			--to-tsv

#produce representative sequence fasta file
qiime feature-table tabulate-seqs  			--i-data rep-seqs.qza  			--o-visualization rep-seqs.qzv
qiime tools export --input-path rep-seqs.qzv  			--output-path unfiltered

#convert to relative abundances
qiime feature-table relative-frequency 			--i-table table.qza 			--o-relative-frequency-table relative-table-ASV.qza

#export to biom
qiime tools export --input-path relative-table-ASV.qza 			--output-path rel-table

#convert to tab separated text file
biom convert 			-i rel-table/feature-table.biom 			-o table/rel-feature-table.tsv --to-tsv
