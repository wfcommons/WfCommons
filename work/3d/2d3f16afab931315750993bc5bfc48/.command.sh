#!/bin/bash -euo pipefail
head -n 1 manifest.txt > header.txt
tail -n+2 manifest.txt | cut -d, -f1 > col1.txt
tail -n+2 manifest.txt | cut -d, -f2 | sed 's:.*/::' > col2.txt
while read f; do
	realpath $f >> full_path.txt
done <col2.txt
tail -n+2 manifest.txt | cut -d, -f3 > col3.txt
paste -d, col1.txt full_path.txt col3.txt > cols.txt
cat cols.txt >> header.txt && mv header.txt manifest.txt

qiime tools import 				--type 'SampleData[PairedEndSequencesWithQuality]' 				--input-path manifest.txt 				--output-path demux.qza 				--input-format PairedEndFastqManifestPhred33
